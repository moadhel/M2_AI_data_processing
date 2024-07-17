% Estimation d'une matrice d'homographie H qui permet de passer d'une 
% image I1 a une autre image I2 a partir de paires de points (homologues) 
% 
% H = [ h11 h12 h13 ; h21 h22 h23 ; h31 h32 h33 ]
% H possede 8 parametres independants. 
% Chaque correspondance donne 2 equations. 
% Ainsi pour estimer H, il faut au moins 4 paires de points. 
%
% Il existe differentes manieres d'estimer H. 
% Nous choisissons la resolution sous la contrainte ||h3|| = 1, 
% au sens des moindres carres. 

function [H] = homographie(XY_C1,XY_C2)
% Entrees :
%
% XY_C1 : matrice (NbPointsx2) contenant les coordonnees des Nbpoints dans l'image I1
% XY_C2 : matrice (NbPointsx2) contenant les coordonnees des Nbpoints HOMOLOGUES dans l'image I2 
%	(colonne 1 : les x, colonne 2 : les y)
%
% Sortie :
% H : la matrice d'homographie estimee


% Les parametres hij de la matrice d'homographie H sont ranges dans 
% un vecteur : H = [ h11 ... h33 ]' tel que 
%           A * H = 0
% avec A qui depend des coordonnees de paires homologues, cf. equation (2)
% A = ( XY_C1(1,1) XY_C1(1,2) 1 0          0          0 -XY_C1(1,1)*XY_C2(1,1) -XY_C1(1,2)*XY_C2(1,1) -XY_C2(1,1) 
%       0          0          0 XY_C1(1,1) XY_C1(1,2) 1 -XY_C1(1,1)*XY_C2(1,2) -XY_C1(1,2)*XY_C2(1,2) -XY_C2(1,2)
%       ... etc ... )

% Stocker dans une variable le nombre de points apparies
NbPoints=size(XY_C1,1);
% Construction des matrices/vecteurs utiles pour construire la matrice A

constr1=zeros(2*NbPoints,1);
constr1(1:2:end)=XY_C1(:,1);

constr2=zeros(2*NbPoints,1);
constr2(1:2:end)=XY_C1(:,2);

constr3=zeros(2*NbPoints,1);
constr3(1:2:end)=ones(NbPoints,1);

constr4=zeros(2*NbPoints,1);
constr4(2:2:end)=XY_C1(:,1);

constr5=zeros(2*NbPoints,1);
constr5(2:2:end)=XY_C1(:,2);

constr6=zeros(2*NbPoints,1);
constr6(2:2:end)=ones(NbPoints,1);

constr7=zeros(2*NbPoints,1);
constr7(1:2:end)=-XY_C1(:,1).*XY_C2(:,1);

constr7(2:2:end)=-XY_C1(:,1).*XY_C2(:,2);

constr8=zeros(2*NbPoints,1);
constr8(1:2:end)=-XY_C1(:,2).*XY_C2(:,1);

constr8(2:2:end)=-XY_C1(:,2).*XY_C2(:,2);

constr9=zeros(2*NbPoints,1);
constr9(1:2:end)=-XY_C2(:,1);

constr9(2:2:end)=-XY_C2(:,2);

A=[constr1 constr2 constr3 constr4 constr5 constr6 constr7 constr8 constr9];

% Estimation des parametres de H par decomposition en valeurs singulieres
% Utiliser la fonction matlab svd : 
% H est le vecteur propre associee a la plus petite valeur propre de A^TA
[U,S,~]=svd(A'*A);

[~, ind]=min(S(S>0),[],1);

H=U(:,ind);




% Former la matrice H de taille 3x3
H=[H(1) H(2) H(3);H(4) H(5) H(6); H(7) H(8) H(9)];