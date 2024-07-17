load("Bmode.mat")
dicomread("20240108163412___.dcm");
figure()
colormap("gray")
imagesc(log(5+6*double(ans(106:1384,604:1473,1))))
figure()
colormap("gray")
imagesc(imgaussfilt(log(5+6*Bmode_invitro),2))
[nRows, nCols] = size(Bmode_invitro);
I_next = log(5+6*double(ans(106:1384,604:1473,1))); 
lambda=0.10;
k = 35; 
num_iterations =30; 
for t = 1:num_iterations
    [Ix, Iy] = gradient(I_next);

    
    cN = 1 ./ (1 + Ix.^2 / k^2); 
    cS = 1 ./ (1 + Iy.^2 / k^2); 

    
    div_cI = cN.*Ix + cS.*Iy; 
    
    % Mise à jour de l'image pour le prochain pas de temps
    I_next = I_next + lambda*div_cI;
end

% Afficher l'image résultante
figure();
imagesc(I_next);
colormap(gray);

%%
hist(exp(I_next(:)),64)
    