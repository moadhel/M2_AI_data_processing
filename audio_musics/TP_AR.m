%% Partie 1
% 1)
% Étape 1 : Charger le fichier audio
file_path = 'segment.dat';
signal = load(file_path);
[signal1,fs]=audioread('signal.wav');
% Visualiser le signal
figure;
plot(signal1);
title('Signal audio');
xlabel('Échantillons');
ylabel('Amplitude');
%% Etape 2 : Calcul du spectre de puissance
L = length(signal1);
S = abs(fft(signal1)).^2 / L;
S = S(1:L/2+1);
% Afficher le logarithme du spectre de puissance
figure;
plot(10*log10(S));
title('Logarithme du spectre de puissance (périodogramme)');
xlabel('Fréquence');
ylabel('Puissance (dB)');
% 2)
%% Etape 1 : Spécifiez l'ordre du modèle AR
P = 200;
% Utilisez la fonction aryule pour estimer les coefficients AR
[a, v_e] = aryule(signal, P);
a = a.';
% Afficher les coefficients AR
disp('Coefficients AR estimés :');
disp(a);
%% Etape 2 : Calcul de la DSP du modèle AR
D = v_e./abs(fft(a,L)).^2;
D = D(1:L/2+1);
% Afficher le logarithme de la DSP
figure;
plot(10*log10(D));
title('Logarithme de la Densité Spectrale de Puissance (DSP) du modèle AR');
xlabel('Fréquence');
ylabel('Puissance (dB)');
%% Etape 3 :
% Fréquences
frequencies = (0:length(signal)/2) / length(signal);
% Affichage des logarithmes du périodogramme et de la DSP
figure;
subplot(2, 1, 1);
plot(frequencies, 10*log10(S));
title('Logarithme du Périodogramme');
xlabel('Fréquence');
ylabel('Puissance (dB)');
subplot(2, 1, 2);
plot(frequencies, 10*log10(D));
title('Logarithme de la DSP du Modèle AR');
xlabel('Fréquence');
ylabel('Puissance (dB)');
sgtitle('Comparaison entre Périodogramme et DSP du Modèle AR');
% Erreur de prédiction
err = filter(a,1,signal) ;
% Afficher l'erreur de prédiction
figure;
subplot(3, 1, 1);
plot(err);
title('Erreur de Prédiction (Résidus)');
xlabel('Échantillons');
ylabel('Amplitude');
% Calculer la variance de l'erreur de prédiction
variance_erreur = var(err);
disp(['Variance de l''erreur de prédiction : ', num2str(variance_erreur)]);
% Calculer et afficher l'auto-corrélation de l'erreur de prédiction
[auto_corr, lags] = xcorr(err, 'coeff');
subplot(3, 1, 2);
plot(lags, auto_corr);
title('Auto-corrélation de l''Erreur de Prédiction');
xlabel('Lags');
ylabel('Auto-corrélation');
% Afficher la variance théorique des résidus (erreur de prédiction)
theoretical_variance = v_e; % e est la variance des résidus de la fonction aryule
disp(['Variance théorique des résidus (erreur de prédiction) : ', num2str(theoretical_variance)]);
sgtitle('Analyse de l''Erreur de Prédiction');
%% Partie 2
% Paramètres du défaut
t01 = 1001;  % Échantillon de départ du défaut
t02=1001;
L_al = 1000;  % Durée du défaut en nombre d'échantillons
amp = 0.5;
% Créer un signal altéré avec le défaut
x = signal1 ;
% Afficher le signal original et le signal altéré
figure;
subplot(2, 1, 1);
plot(signal);
title('Signal Original');
xlabel('Échantillons');
ylabel('Amplitude');
subplot(2, 1, 2);
plot(x);
title('Signal Altéré avec Défaut');
xlabel('Échantillons');
ylabel('Amplitude');
% Ordre du filtre AR
P = 40;  % Remplacez par l'ordre souhaité
% Estimer les coefficients AR et la variance des résidus
[a, v_e] = aryule(x, P);
a = a.';
% Filtrer le signal altéré pour obtenir l'erreur de prédiction
erreur_prediction = filter(a,1, x);
% Afficher le signal original, le signal altéré et l'erreur de prédiction
figure;
subplot(3, 1, 1);
plot(signal);
title('Signal Original');
xlabel('Échantillons');
ylabel('Amplitude');
subplot(3, 1, 2);
plot(x);
title('Signal Altéré avec Défaut');
xlabel('Échantillons');
ylabel('Amplitude');
subplot(3, 1, 3);
plot(erreur_prediction);
title('Erreur de Prédiction (Résidus)');
xlabel('Échantillons');
ylabel('Amplitude');
sgtitle('Analyse du Signal Altéré avec Modèle AR');
% Afficher la variance des résidus
disp(['Variance des résidus (erreur de prédiction) : ', num2str(v_e)]);
% Calculer et afficher l'auto-corrélation de l'erreur de prédiction
[auto_corr, lags] = xcorr(erreur_prediction, 'coeff');
figure;
plot(lags, auto_corr);
title('Auto-corrélation de l''Erreur de Prédiction');
xlabel('Lags');
ylabel('Auto-corrélation');
% détecteur de clics
detecteur = erreur_prediction > 3*sqrt(v_e);
% Afficher le signal altéré et le détecteur
figure;
subplot(2, 1, 1);
plot(x);
title('Signal Altéré avec Défaut');
xlabel('Échantillons');
ylabel('Amplitude');
subplot(2, 1, 2);
plot(detecteur, 'r');
title('Détecteur (Rouge : Détection Anomalie)');
xlabel('Échantillons');
ylabel('Détection');
sgtitle('Détection d''Anomalie avec Détecteur');
%% étape 1 2 3
P=500;
L=3000;
s_al1=zeros(1,L_al);
s_al2=zeros(1,L_al);

for i=1:3
    [a, v_e] = aryule(x(27000:30000), P);
    R = zeros(1,L);
    R(1:P+1) = a(end:-1:1);
    C = zeros(L-P,1);
    C(1) = a(P+1);
    A = toeplitz(C,R);
    % etape 4
    Aav=A(:,1:t01-1);
    Aal=A(:,t01:t01+L_al-1);
    Aap=A(:,t01+L_al:end);
    % étape 5
    s_al1=-(Aal'*Aal)\Aal'*(Aav*x(27001:28000)+Aap*x(29001:30000));

end
for i=1:3
    [a, v_e] = aryule(x(52000:55000), P);
    R = zeros(1,L);
    R(1:P+1) = a(end:-1:1);
    C = zeros(L-P,1);
    C(1) = a(P+1);
    A = toeplitz(C,R);
    % etape 4
    Aav=A(:,1:t02-1);
    Aal=A(:,t02:t02+L_al-1);
    Aap=A(:,t02+L_al:end);
    % étape 5
    s_al2=-(Aal'*Aal)\Aal'*(Aav*x(52001:53000)+Aap*x(54001:55000));

end
%% partie 3
s_hat=[x(1:27999)' s_al1' x(29001:52999)' s_al2' x(54001:end)'];

figure;
subplot(2, 1, 1);
plot(x);
title('Signal Altéré avec Défaut');
xlabel('Échantillons');
ylabel('Amplitude');
subplot(2, 1, 2);
plot(s_hat, 'r');
title('Signal reconstruit');
audiowrite('signal_ss_clique.wav', s_hat,44100);

