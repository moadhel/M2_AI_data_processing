function diff_im = anisodiff2D(im, num_iter, lambda, k, option)
% Convert input image to double.
im = double(im);

% Condition initiale
diff_im = im;

% Masques pour le calcul des gradients.
hW = [0 0 0; 1 -1 0; 0 0 0];
hS = rot90(hW);
hE = rot90(hS);
hN = rot90(hE);

% Anisotropic diffusion.
for t = 1:num_iter

% Calcul des gradients par convolution.
nablaN = imfilter(diff_im,hN,'conv');
nablaS = imfilter(diff_im,hS,'conv');   
nablaW = imfilter(diff_im,hW,'conv');
nablaE = imfilter(diff_im,hE,'conv');   

% Expression du coefficient de diffusion.
if option == 1
cN = exp(-(nablaN/k).^2);
cS = exp(-(nablaS/k).^2);
cW = exp(-(nablaW/k).^2);
cE = exp(-(nablaE/k).^2);

elseif option == 2
cN = 1./(1 + (nablaN/k).^2);
cS = 1./(1 + (nablaS/k).^2);
cW = 1./(1 + (nablaW/k).^2);
cE = 1./(1 + (nablaE/k).^2);

end

% Mise à jour de la solution itérative
diff_im = diff_im + lambda*(cN.*nablaN+cS.*nablaS+cW.*nablaW+cE.*nablaE);

end

end