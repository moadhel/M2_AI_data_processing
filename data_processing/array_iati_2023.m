%Array processing course basic code
clear
clc
close all
format shortG
%+++++ BEAMFORMING ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%----- Scenario -----
%Number of elements in the array
N = 10;
%Inter-element spacing (in wavelength)
d = 0.5;
pos = d * (0:N-1)'; %positions of the antennas
%Mainlobe width
theta_3dB = 0.9/(N*d);
%White noise
sigma2 = 1;	%white noise power
%Interference
thetaj = [-20;15]/180*pi;	%angles of arrival	
INR = [20;20];			%interference to noise ratio (dB)
Pj = sigma2 * 10.^(INR/10);		%interference power
J = length(thetaj);
%Interference + noise covariance matrix
Aj = exp(1i*2*pi*pos*sin(thetaj'));	%interference steering matrix N|J
C = Aj*diag(Pj)*Aj' + sigma2*eye(N);	%interference + noise covariance matrix
%Signal of interest
thetas = 0/180*pi;	%angle of arrival
%thetas = 30/180*pi; %lobe principal décalé à 30°
SNR = 0;            %signal to noise ratio (dB)
Ps = sigma2 * 10^(SNR/10);			%signal power
as = exp(1i*2*pi*pos*sin(thetas));	%steering vector
%Total covariance matrix (signal + interference + noise)
R = Ps*(as*as') + C;


%----- Natural beampattern -----
% %Weight vector (all weights equal to 1/N)
w = ones(N,1); 
w = w/(ones(1,N)*w); 
% w = as; %pour viser une direction suivant thetas
% w = w/(as'*w); %pour normaliser 
%Diagram
tab_theta = (-90:0.5:90)/180*pi;        %Angles where to evaluate beampattern
A = exp(1i*2*pi*pos*sin(tab_theta));    %Steering matrix: each column is a(theta)
G = 20*log10(abs(w'*A));        %beampattern (power in dB)
%Plot
figure
plot(tab_theta*180/pi,G,'linewidth',2);
title('Natural beampattern $w_n=1/N$','interpreter','latex');
ylabel('dB','interpreter','latex');
xlabel('Angle of Arrival (degrees)','interpreter','latex');
axis([-90 90 -70 10]);
grid on


%----- CONVENTIONAL AND OPTIMAL BEAMFORMERS -----
%Looked direction
theta0 = 0/180*pi;
a0 = exp(1i*2*pi*pos*sin(theta0));
%Conventional beamformer
w_CBF = a0; 
w_CBF = w_CBF/(a0'*w_CBF);
G_CBF = 20*log10(abs(w_CBF'*A));        %beampattern (power in dB)
SINR_CBF = Ps*(abs(w_CBF'*as)^2)/(abs(w_CBF'*C*w_CBF)); %SINR
A_WN_CBF = 1/(norm(w_CBF)^2);   %White noise array gain
%Optimal beamformer
w_opt = (C\as); 
w_opt = w_opt/(as'*w_opt);
G_opt = 20*log10(abs(w_opt'*A));
SINR_opt = Ps*(abs(w_opt'*as)^2)/(abs(w_opt'*C*w_opt));
A_WN_opt = 1/(norm(w_opt)^2);
%Plot
figure
plot(tab_theta*180/pi,G_opt,'-',...
    tab_theta*180/pi,G_CBF,'--','linewidth',2);
for k=1:length(thetaj)
    xline(thetaj(k)*180/pi,'--','color','k','linewidth',2)
end
title('Beampatterns of CBF and optimal beamformers','interpreter','latex');
ylabel('dB','interpreter','latex');
xlabel('Angle of Arrival (degrees)','interpreter','latex');
legend('opt','CBF');
axis([-90 90 -70 10]);
grid on


%%
%----- ADAPTIVE BEAMFORMING WITH TRUE COVARIANCE MATRICES -----
%MVDR
w_MVDR = (C\a0); 
w_MVDR = w_MVDR/(a0'*w_MVDR);
G_MVDR = 20*log10(abs(w_MVDR'*A));
SINR_MVDR = Ps*(abs(w_MVDR'*as)^2)/(abs(w_MVDR'*C*w_MVDR));
A_WN_MVDR = 1/(norm(w_MVDR)^2);
%MPDR
w_MPDR = (R\a0); 
w_MPDR = w_MPDR/(a0'*w_MPDR);
G_MPDR = 20*log10(abs(w_MPDR'*A));
SINR_MPDR = Ps*(abs(w_MPDR'*as)^2)/(abs(w_MPDR'*C*w_MPDR));
A_WN_MPDR = 1/(norm(w_MPDR)^2);
%Plot
figure
plot(tab_theta*180/pi,G_opt,'-',...
    tab_theta*180/pi,G_CBF,'--',...
    tab_theta*180/pi,G_MVDR,'-.',...
    tab_theta*180/pi,G_MPDR,':.','linewidth',2);
for k=1:length(thetaj)
    xline(thetaj(k)*180/pi,'--','color','k','linewidth',2)
end
title('Beampatterns $K=\infty$','interpreter','latex');
ylabel('dB','interpreter','latex');
xlabel('Angle of Arrival (degrees)','interpreter','latex');
legend('opt','CBF','MVDR','MPDR');
axis([-90 90 -70 10]);
grid on


disp(['SINR_CBF = ',num2str(10*log10(SINR_CBF)),' dB'])
disp(['SINR_opt = ',num2str(10*log10(SINR_opt)),' dB'])
disp(['SINR_MVDR = ',num2str(10*log10(SINR_MVDR)),' dB'])
disp(['SINR_MPDR = ',num2str(10*log10(SINR_MPDR)),' dB'])

% %----- ADAPTIVE BEAMFORMING WITH ESTIMATED COVARIANCE MATRICES -----
% %Number of snapshots
% 
% T= 1000;
% 
% A_WN_MPDR_1=zeros(1,T);
% A_WN_MVDR_l=zeros(1,T);
% SINR_MPDR_1=zeros(1,T);
% SINR_MVDR_l=zeros(1,T);
% 
% A_WN_MPDR_2=zeros(1,T);
% SINR_MPDR_2=zeros(1,T);
% 
% 
% A_WN_MPDR_3=zeros(1,T);
% SINR_MPDR_3=zeros(1,T);
% 
% 
% mu_tout= [12.4,25,200] ;
% 
% 
%     for K = 1:T
%         %Signal
%         S = sqrt(Ps/2) * as * (randn(1,K)+1i*randn(1,K));
%         %Interference + noise
%         IN = Aj * diag(sqrt(Pj/2)) * (randn(J,K)+1i*randn(J,K));
%         NOISE = sqrt(sigma2/2)*(randn(N,K)+1i*randn(N,K));
%         %MVDR-SMI
%         Y_MVDR = IN + NOISE;
%         C_hat = (Y_MVDR*Y_MVDR')/K;
%         w_MVDR_SMI = (C_hat\a0);
%         w_MVDR_SMI = w_MVDR_SMI / (a0'*w_MVDR_SMI);
%         G_MVDR_SMI = 20*log10(abs(w_MVDR_SMI'*A));
%         SINR_MVDR_SMI = Ps*(abs(w_MVDR_SMI'*as)^2)/(abs(w_MVDR_SMI'*C*w_MVDR_SMI));
%         A_WN_MVDR_SMI = 1 / (norm(w_MVDR_SMI)^2);
%         A_WN_MVDR_l(K)=A_WN_MVDR_SMI;
%         SINR_MVDR_l(K)=SINR_MVDR_SMI;
%         %MPDR-SMI    SINR_MVDR_l(K)=SINR_MVDR_SMI;
%         Y_MPDR = S + IN + NOISE;
%         R_hat = (Y_MPDR*Y_MPDR')/K;
%         %w_MPDR_SMI = (R_hat\a0);
%         w_MPDR_SMI1 = ((R_hat+mu_tout(1)*eye(size(R_hat)))\a0); 
%         w_MPDR_SMI1 = w_MPDR_SMI1 / (a0'*w_MPDR_SMI1);
%         G_MPDR_SMI1 = 20*log10(abs(w_MPDR_SMI1'*A));
%         SINR_MPDR_SMI1 = Ps*(abs(w_MPDR_SMI1'*as)^2)/(abs(w_MPDR_SMI1'*C*w_MPDR_SMI1));
%         A_WN_MPDR_SMI1 = 1 / (norm(w_MPDR_SMI1)^2);
% 
%         w_MPDR_SMI2 = ((R_hat+mu_tout(2)*eye(size(R_hat)))\a0); 
%         w_MPDR_SMI2 = w_MPDR_SMI2 / (a0'*w_MPDR_SMI2);
%         G_MPDR_SMI2 = 20*log10(abs(w_MPDR_SMI2'*A));
%         SINR_MPDR_SMI2 = Ps*(abs(w_MPDR_SMI2'*as)^2)/(abs(w_MPDR_SMI2'*C*w_MPDR_SMI2));
%         A_WN_MPDR_SMI2 = 1 / (norm(w_MPDR_SMI2)^2);
% 
%         w_MPDR_SMI3 = ((R_hat+mu_tout(3)*eye(size(R_hat)))\a0); 
%         w_MPDR_SMI3 = w_MPDR_SMI3 / (a0'*w_MPDR_SMI3);
%         G_MPDR_SMI3 = 20*log10(abs(w_MPDR_SMI3'*A));
%         SINR_MPDR_SMI3 = Ps*(abs(w_MPDR_SMI3'*as)^2)/(abs(w_MPDR_SMI3'*C*w_MPDR_SMI3));
%         A_WN_MPDR_SMI3 = 1 / (norm(w_MPDR_SMI3)^2);
% 
%         A_WN_CBF = 1 / (norm(w_CBF)^2);
%         A_WN_opt = 1 / (norm(w_opt)^2);
%         A_WN_MPDR_1(K)=A_WN_MPDR_SMI1;
%         SINR_MPDR_1(K)=SINR_MPDR_SMI1;
%         A_WN_MPDR_2(K)=A_WN_MPDR_SMI2;
%         SINR_MPDR_2(K)=SINR_MPDR_SMI2;
%         A_WN_MPDR_3(K)=A_WN_MPDR_SMI3;
%         SINR_MPDR_3(K)=SINR_MPDR_SMI3;
%     end
% 
% p = 50;
% 
% for i = 1:(T-p)
%     moy1 = sum(A_WN_MPDR_1(i:i+p))/(p+1);
%     moy2= sum(A_WN_MVDR_l(i:i+p))/(p+1);
%     moy3 = sum(SINR_MPDR_1(i:i+p))/(p+1);
%     moy4 = sum(SINR_MVDR_l(i:i+p))/(p+1);
%     moy5= sum(A_WN_MPDR_2(i:i+p))/(p+1);
%     moy6 = sum(SINR_MPDR_2(i:i+p))/(p+1);
%     moy7= sum(A_WN_MPDR_3(i:i+p))/(p+1);
%     moy8 = sum(SINR_MPDR_3(i:i+p))/(p+1);
%     for l = 1:p
%         A_WN_MPDR_1(i+l)=moy1;
%         A_WN_MVDR_l(i+l)=moy2;
%         SINR_MPDR_1(i+l)=moy3;
%         SINR_MVDR_l(i+l)=moy4;
%         A_WN_MPDR_2(i+l)=moy5;
%         SINR_MPDR_2(i+l)=moy6;
%         A_WN_MPDR_3(i+l)=moy7;
%         SINR_MPDR_3(i+l)=moy8;
%     end
% end
% 
% % %Plot
% % figure
% % plot(tab_theta*180/pi,G_opt,'-',...
% %     tab_theta*180/pi,G_CBF,'--',...
% %     tab_theta*180/pi,G_MVDR_SMI,'-.',...
% %     tab_theta*180/pi,G_MPDR_SMI,':.','linewidth',2);
% % for k=1:length(thetaj)
% %     xline(thetaj(k)*180/pi,'--','color','k','linewidth',2)
% % end
% % title(['Beampatterns $K=$',num2str(K)],'interpreter','latex');
% % ylabel('dB','interpreter','latex');
% % xlabel('Angle of Arrival (degrees)','interpreter','latex');
% % legend('opt','CBF','MVDR-SMI','MPDR-SMI');
% % axis([-90 90 -70 10]);
% % grid on
% % 
% % 
% % disp(['SINR_MVDR_SMI = ',num2str(10*log10(SINR_MVDR_SMI)),' dB'])
% % disp(['SINR_MPDR_SMI = ',num2str(10*log10(SINR_MPDR_SMI)),' dB'])
% % %White noise array gain
% % disp(['A_WN_CBF = ',num2str(10*log10(A_WN_CBF)),' dB'])
% % disp(['A_WN_opt = ',num2str(10*log10(A_WN_opt)),' dB'])
% % disp(['A_WN_MVDR = ',num2str(10*log10(A_WN_MVDR)),' dB'])
% % disp(['A_WN_MPDR = ',num2str(10*log10(A_WN_MPDR)),' dB'])
% % disp(['A_WN_MVDR_SMI = ',num2str(10*log10(A_WN_MVDR_SMI)),' dB'])
% % disp(['A_WN_MPDR_SMI = ',num2str(10*log10(A_WN_MPDR_SMI)),' dB'])
% 
% figure
% hold on
% plot(1:T,A_WN_MPDR_1)
% plot(1:T,SINR_MPDR_1)
% plot(1:T,A_WN_MPDR_2)
% plot(1:T,SINR_MPDR_2)
% plot(1:T,A_WN_MPDR_3)
% plot(1:T,SINR_MPDR_3)
% % plot(1:T,A_WN_MVDR_l)
% % plot(1:T,SINR_MVDR_l)
% legend('A_W_N MPDR1','SINR MPDR1','A_W_N MPDR2','SINR MPDR2','A_W_N MPDR3','SINR MPDR3')
% hold off

% %----- ADAPTIVE BEAMFORMING WITH ESTIMATED COVARIANCE MATRICES -----
% %Number of snapshots
% 
% T= 1000;
% 
% A_WN_MPDR_l=zeros(1,T);
% A_WN_MVDR_l=zeros(1,T);
% SINR_MPDR_l=zeros(1,T);
% SINR_MVDR_l=zeros(1,T);
% 
% mu = 0;
% 
% 
% for K = 1:T
%     %Signal
%     S = sqrt(Ps/2) * as * (randn(1,K)+1i*randn(1,K));
%     %Interference + noise
%     IN = Aj * diag(sqrt(Pj/2)) * (randn(J,K)+1i*randn(J,K));
%     NOISE = sqrt(sigma2/2)*(randn(N,K)+1i*randn(N,K));
%     %MVDR-SMI
%     Y_MVDR = IN + NOISE;
%     C_hat = (Y_MVDR*Y_MVDR')/K;
%     w_MVDR_SMI = (C_hat\a0);
%     w_MVDR_SMI = w_MVDR_SMI / (a0'*w_MVDR_SMI);
%     G_MVDR_SMI = 20*log10(abs(w_MVDR_SMI'*A));
%     SINR_MVDR_SMI = Ps*(abs(w_MVDR_SMI'*as)^2)/(abs(w_MVDR_SMI'*C*w_MVDR_SMI));
%     A_WN_MVDR_SMI = 1 / (norm(w_MVDR_SMI)^2);
%     A_WN_MVDR_l(K)=A_WN_MVDR_SMI;
%     SINR_MVDR_l(K)=SINR_MVDR_SMI;
%     %MPDR-SMI    SINR_MVDR_l(K)=SINR_MVDR_SMI;
%     Y_MPDR = S + IN + NOISE;
%     R_hat = (Y_MPDR*Y_MPDR')/K;
%     w_MPDR_SMI = (R_hat\a0);
%     %w_MPDR_SMI = ((R_hat+mu*eye(size(R_hat)))\a0); 
%     w_MPDR_SMI = w_MPDR_SMI / (a0'*w_MPDR_SMI);
%     G_MPDR_SMI = 20*log10(abs(w_MPDR_SMI'*A));
%     SINR_MPDR_SMI = Ps*(abs(w_MPDR_SMI'*as)^2)/(abs(w_MPDR_SMI'*C*w_MPDR_SMI));
%     A_WN_MPDR_SMI = 1 / (norm(w_MPDR_SMI)^2);
%     A_WN_CBF = 1 / (norm(w_CBF)^2);
%     A_WN_opt = 1 / (norm(w_opt)^2);
%     A_WN_MPDR_l(K)=A_WN_MPDR_SMI;
%     SINR_MPDR_l(K)=SINR_MPDR_SMI;
% end
% 
% p = 100;
% 
% for i = 1:(T-p)
%     moy1 = sum(A_WN_MPDR_l(i:i+p))/(p+1);
%     moy2= sum(A_WN_MVDR_l(i:i+p))/(p+1);
%     moy3 = sum(SINR_MPDR_l(i:i+p))/(p+1);
%     moy4 = sum(SINR_MVDR_l(i:i+p))/(p+1);
%     for l = 1:p
%         A_WN_MPDR_l(i+l)=moy1;
%         A_WN_MVDR_l(i+l)=moy2;
%         SINR_MPDR_l(i+l)=moy3;
%         SINR_MVDR_l(i+l)=moy4;
%     end
% end
% 
% %Plot
% figure
% plot(tab_theta*180/pi,G_opt,'-',...
%     tab_theta*180/pi,G_CBF,'--',...
%     tab_theta*180/pi,G_MVDR_SMI,'-.',...
%     tab_theta*180/pi,G_MPDR_SMI,':.','linewidth',2);
% for k=1:length(thetaj)
%     xline(thetaj(k)*180/pi,'--','color','k','linewidth',2)
% end
% title(['Beampatterns $K=$',num2str(K)],'interpreter','latex');
% ylabel('dB','interpreter','latex');
% xlabel('Angle of Arrival (degrees)','interpreter','latex');
% legend('opt','CBF','MVDR-SMI','MPDR-SMI');
% axis([-90 90 -70 10]);
% grid on
% 
% 
% disp(['SINR_MVDR_SMI = ',num2str(10*log10(SINR_MVDR_SMI)),' dB'])
% disp(['SINR_MPDR_SMI = ',num2str(10*log10(SINR_MPDR_SMI)),' dB'])
% %White noise array gain
% disp(['A_WN_CBF = ',num2str(10*log10(A_WN_CBF)),' dB'])
% disp(['A_WN_opt = ',num2str(10*log10(A_WN_opt)),' dB'])
% disp(['A_WN_MVDR = ',num2str(10*log10(A_WN_MVDR)),' dB'])
% disp(['A_WN_MPDR = ',num2str(10*log10(A_WN_MPDR)),' dB'])
% disp(['A_WN_MVDR_SMI = ',num2str(10*log10(A_WN_MVDR_SMI)),' dB'])
% disp(['A_WN_MPDR_SMI = ',num2str(10*log10(A_WN_MPDR_SMI)),' dB'])
% 
% figure
% hold on
% plot(1:T,A_WN_MPDR_l, 'r')
% plot(1:T,SINR_MPDR_l, 'g')
% plot(1:T,A_WN_MVDR_l, 'blue')
% plot(1:T,SINR_MVDR_l,'black')
% legend('A_W_N MPDR','SINR MPDR','A_W_N MVDR','SINR MVDR')
% hold off


% 
%----- ADAPTIVE BEAMFORMING WITH ESTIMATED COVARIANCE MATRICES -----
fin=100;
awn = zeros(1,991);
sinrr = zeros(1,991);

mu = 1:0.1:fin;
for y = 1:991
    %Number of snapshots
    K = 30;
    %Signal
    S = sqrt(Ps/2) * as * (randn(1,K)+1i*randn(1,K));
    %Interference + noise
    IN = Aj * diag(sqrt(Pj/2)) * (randn(J,K)+1i*randn(J,K));
    NOISE = sqrt(sigma2/2)*(randn(N,K)+1i*randn(N,K));
    %MVDR-SMI
    Y_MVDR = IN + NOISE;
    C_hat = (Y_MVDR*Y_MVDR')/K;
    w_MVDR_SMI = (C_hat\a0);
    w_MVDR_SMI = w_MVDR_SMI / (a0'*w_MVDR_SMI);
    G_MVDR_SMI = 20*log10(abs(w_MVDR_SMI'*A));
    SINR_MVDR_SMI = Ps*(abs(w_MVDR_SMI'*as)^2)/(abs(w_MVDR_SMI'*C*w_MVDR_SMI));
    A_WN_MVDR_SMI = 1 / (norm(w_MVDR_SMI)^2);
    %MPDR-SMI
    Y_MPDR = S + IN + NOISE;
    R_hat = (Y_MPDR*Y_MPDR')/K;
    w_MPDR_SMI = ((R_hat+mu(y)*eye(size(R_hat)))\a0);
    w_MPDR_SMI = w_MPDR_SMI / (a0'*w_MPDR_SMI);
    G_MPDR_SMI = 20*log10(abs(w_MPDR_SMI'*A));
    SINR_MPDR_SMI = Ps*(abs(w_MPDR_SMI'*as)^2)/(abs(w_MPDR_SMI'*C*w_MPDR_SMI));
    A_WN_MPDR_SMI = 1 / (norm(w_MPDR_SMI)^2);
    A_WN_CBF = 1 / (norm(w_CBF)^2);
    A_WN_opt = 1 / (norm(w_opt)^2);
    
    awn(round(y))=A_WN_MPDR_SMI;
    sinrr(round(y)) = SINR_MPDR_SMI;
end
t=991;
p=100;
for i = 1:(t-p)
    moy1 = sum(awn(i:i+p))/(p+1);
    moy2 = sum(sinrr(i:i+p))/(p+1);
    for l = 1:p
        awn(i+l)=moy1;
        sinrr(i+l)=moy2;
    end
end

%Plot
figure
plot(tab_theta*180/pi,G_opt,'-',...
    tab_theta*180/pi,G_CBF,'--',...
    tab_theta*180/pi,G_MVDR_SMI,'-.',...
    tab_theta*180/pi,G_MPDR_SMI,':.','linewidth',2);
for k=1:length(thetaj)
    xline(thetaj(k)*180/pi,'--','color','k','linewidth',2)
end
title(['Beampatterns $K=$',num2str(K)],'interpreter','latex');
ylabel('dB','interpreter','latex');
xlabel('Angle of Arrival (degrees)','interpreter','latex');
legend('opt','CBF','MVDR-SMI','MPDR-SMI');
axis([-90 90 -70 10]);
grid on


disp(['SINR_MVDR_SMI = ',num2str(10*log10(SINR_MVDR_SMI)),' dB'])
disp(['SINR_MPDR_SMI = ',num2str(10*log10(SINR_MPDR_SMI)),' dB'])
%White noise array gain
disp(['A_WN_CBF = ',num2str(10*log10(A_WN_CBF)),' dB'])
disp(['A_WN_opt = ',num2str(10*log10(A_WN_opt)),' dB'])
disp(['A_WN_MVDR = ',num2str(10*log10(A_WN_MVDR)),' dB'])
disp(['A_WN_MPDR = ',num2str(10*log10(A_WN_MPDR)),' dB'])
disp(['A_WN_MVDR_SMI = ',num2str(10*log10(A_WN_MVDR_SMI)),' dB'])
disp(['A_WN_MPDR_SMI = ',num2str(10*log10(A_WN_MPDR_SMI)),' dB'])

figure 
hold on
plot(1:0.1:fin,awn)
plot(1:0.1:fin,sinrr)
legend('awn','sinrr')
%%
%----- GSC implementation of MVDR beamformer -----
%Matrix B
B = null(a0');
%MVDR
Y = Y_MVDR;
%Data in the main and auxilliary channels
d = w_CBF' * Y;     %signal in main channel 1|K
Z = B' * Y;         %signal in auxilliary channels N-1|K
%Estimate covariance matrix and cross correlation
Rz = (Z*Z')/K;      %estimate of R_z
rdz = Z*d'/K;       %estimate of R_{dz}
wa = Rz\rdz;        %estimate of R_z^{-1} r_{dz}
w_MVDR_SMI_GSC = w_CBF - B * wa;
%disp(['||w_{gsc}-w_{df}||=',num2str(norm(w_MVDR_SMI-w_MVDR_SMI_GSC))])

%----- GSC implementation of MVDR beamformer -----
%Matrix B
[U,S,V]=svd(Rz);
Ps1 = sigma2 * 10^(0/10);

K=1000;
N=10;
w_MVDR_SMI_GSC2=zeros(1,K);
B = null(a0');
%U=B'*Aj;
%MVDR
for R=1:9
%Ps1 = sigma2 * 10^(0/10);
% for i=1:K
% IN = Aj * diag(sqrt(Pj/2)) * (randn(J,i)+1i*randn(J,i));
% NOISE = sqrt(sigma2/2)*(randn(N,i)+1i*randn(N,i));
%MVDR-SMI
% Y_MVDR = IN + NOISE;
Y = Y_MVDR;
%Data in the main and auxilliary channels
d = w_CBF' * Y;     %signal in main channel 1|K
Z = B' * Y; %signal in auxilliary channels N-1|K
Z2=U(:,1:R)'*Z(:,:);

%Estimate covariance matrix and cross correlation
Rz2 = (Z2*Z2')/K;      %estimate of R_z
rdz2 = Z2*d'/K;       %estimate of R_{dz}
wa2 = Rz2\rdz2;        %estimate of R_z^{-1} r_{dz}
w_MVDR_SMI_GSC22 = w_CBF - B* U(:,1:R)* wa2;
% w_MVDR_SMI_GSC2(i)=norm(w_MVDR_SMI-w_MVDR_SMI_GSC22);
disp(['||w_{gsc}-w_{df}||=',num2str(Ps1*(abs(w_MVDR_SMI_GSC22'*as)^2)/(abs(w_MVDR_SMI_GSC22'*C*w_MVDR_SMI_GSC22)))])
%disp(Ps1*(abs(w_MVDR_SMI_GSC22'*as)^2)/(abs(w_MVDR_SMI_GSC22'*C*w_MVDR_SMI_GSC22)))
end
%+++++++++++ DIRECTION OF ARRIVAL ESTIMATION ++++++++++++++++++++++++++++++
%Number of elements in the array
N = 20;
%Inter-element spacing (in wavelength)
d = 0.5;
pos = d * (0:N-1)'; %positions of the antennas
%Signals impinging on the array 
%White noise
sigma2 = 1;	%white noise power
%Signals
thetas = [-20;0;10]/180*pi;	%angles of arrival	
SNR = [10;10;10];			%interference to noise ratio
Ps = sigma2 * 10.^(SNR/10);		%interference power
P = length(thetas);
As = exp(1i*2*pi*pos*sin(thetas'));
%Number of snapshots
K = N;
%Number of FFT bins
vec_doa = (-90:0.2:90)/180*pi;       
A = exp(1i*2*pi*pos*sin(vec_doa)); 
%Snapshots
Y = As * diag(sqrt(Ps/2)) * (randn(P,K)+1i*randn(P,K))...
    + sqrt(sigma2/2)*(randn(N,K)+1i*randn(N,K));

%Conventional beamforming
P_CBF = sum(abs(A'*Y).^2,2)/(N^2)/K;
P_CBF = 10*log10(P_CBF);

figure
vec_doa = vec_doa*180/pi;
plot(vec_doa,P_CBF,'-','Linewidth',2);
title('CBF spectrum','interpreter','latex');
xlabel('Angle of arrival (degrees)','interpreter','latex');
ylabel('dB','interpreter','latex');

