clear
clc
close all
format shortG
%+++++ SECTIONS 4 AND 5.1 +++++++++++++++++++++++++++++++++++++++++++++++++
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
SNR = 0;            %signal to noise ratio (dB)
Ps = sigma2 * 10^(SNR/10);			%signal power
as = exp(1i*2*pi*pos*sin(thetas));	%steering vector
%Total covariance matrix (signal + interference + noise)
R = Ps*(as*as') + C;

%----- AUXILIARY CHANNELS -------------------------------------------------
%Blocking matrix
a0 = as;
B = null(a0');
tab_theta = (-90:0.5:90)/180*pi;        %Angles where to evaluate beampattern
A = exp(1i*2*pi*pos*sin(tab_theta));    %Steering matrix: each column is a(theta)
G_B = 20*log10(abs(B'*A)).';
figure
plot(tab_theta*180/pi,G_B,'linewidth',2);
for k=1:length(thetaj)
    xline(thetaj(k)*180/pi,'--','color','k','linewidth',2)
end
title('Beampatterns of the columns of $\mathbf{B}$','interpreter','latex');
ylabel('dB','interpreter','latex');
xlabel('Angle of Arrival (degrees)','interpreter','latex');
axis([-90 90 -30 10]);
%Test another matrix B
%Orthogonal (N-1)|(N-1) matrix
U = orth(randn(N-1,N-1) + 1i*randn(N-1,N-1));
B2 = B*U;
G_B2 = 20*log10(abs(B2'*A)).';
figure
plot(tab_theta*180/pi,G_B2,'linewidth',2);
for k=1:length(thetaj)
    xline(thetaj(k)*180/pi,'--','color','k','linewidth',2)
end
title('Beampatterns of the columns of $\mathbf{B}_{2}$','interpreter','latex');
ylabel('dB','interpreter','latex');
xlabel('Angle of Arrival (degrees)','interpreter','latex');
axis([-90 90 -30 10]);


%----- CHECK THAT OPTIMAL BEAMFORMER IS PARTIALLY ADAPTIVE -----
%MVDR beamformer (known C, aO=as)
a0 = as;
w_MVDR = (C\a0); 
w_MVDR = w_MVDR/(a0'*w_MVDR);
%Blocking matrix
B = null(a0');
%Check that w_MVDR belongs to subspace spanned by [a0 , B B^H Aj]
[Q,~] = qr([a0 B*B'*Aj],0); %orthogonal basis Q^H Q = I_{J+1}
disp(norm(Q'*Q-eye(J+1)));
disp(norm(w_MVDR-Q*Q'*w_MVDR));