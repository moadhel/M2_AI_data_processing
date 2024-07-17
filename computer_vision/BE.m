%% Partie I

load('object.mat'); % Lecture objet u
load('psf.mat'); % Lecture psf k

imdisp(u,'Objet');
imdisp(fftshift(k),'PSF');
imdisp(conv2((k),(u)),'Image acquisition')

imdisp(log(1+abs(fftshift(fft2(k).*fft2(u)))),'Image dans le domaine de Fourier')

%% Partie II
load('pattern_single.mat'); % Lecture illumination w
imdisp(w,'Illumination w');
imdisp(log(1+abs(fft2(w))),'Illumination Transformé de Fourier w');
imdisp(w.*u,'Objet super résolu');

imdisp(log(1+abs(fftshift(fft2(k).*fft2(w.*u)))),'Image super résolu')
%% Partie III

load('acquisitions.mat'); % Lecture des trois acquisitions v
load('pattern_full.mat'); % Lecture des trois illuminations w_full
figure;
for ii=1:3
   
    subplot(2,3,ii);
    imdisp(v(:,:,ii),['Acquisition #',num2str(ii)],0);
    subplot(2,3,ii+3);
    imdisp(w_full(:,:,ii),['Pattern #',num2str(ii)],0);
end
alpha=0.01;
lambda =0.01; 
max_iterations = 100; 

[M, N] = size(v(:,:,1));

u_reconstructed = ones(M, N);

for iteration = 1:max_iterations
    gradient = zeros(M, N);
    for ii = 1:3
        v_reconstructed=ifft2((fft2(k).*fft2(w_full(:,:,ii).*u_reconstructed)));
        gradient = gradient - 2 *  (w_full(:,:,ii)) .*ifft2(conj(fft2(k)).*fft2(v(:,:,ii)))-w_full(:,:,ii).*ifft2(abs(fft2(k)).^2.*fft2(w_full(:,:,ii).*u_reconstructed))+alpha*u_reconstructed;
    end
    
    u_reconstructed = u_reconstructed - lambda * gradient;
    
    u_reconstructed(u_reconstructed < 0) = 0;
end

imdisp(u_reconstructed, 'Image reconstruite');
imdisp(log(abs(1+fftshift(fft2(u_reconstructed)))), 'Transformée de Fourier de l''image reconstruite');

