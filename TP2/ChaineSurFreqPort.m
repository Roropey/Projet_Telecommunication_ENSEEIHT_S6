%% Nettoyage
close all;
clear;

%% Variables Initiales
alpha  = 0.35;
Fp = 2000;
Fe = 10000;
Rb = 2000;
nb_bits = 10000;
info_binaire = randi([0,1], 1,nb_bits);

%% Modulateur

% Variables
Ns = (Fe/Rb)*2;
a_00 = 0;
a_01 = 1;
a_11 = 2;
a_10 = 3;
h = ones(1,Ns);
M = 4;
t = [0:1/Fe:1/Fe*nb_bits/2];

% Calculs
info_binaire_2 = reshape(info_binaire, [2 nb_bits/2]);
mapping = info_binaire_2(1,:).*info_binaire_2(2,:).*(a_00-a_01-a_10+a_11) + info_binaire_2(1,:).*(a_01-a_11) + info_binaire_2(2,:).*(a_10-a_11) + a_11;
Suite_diracs = kron(mapping);
xe = filter(h, 1, Suite_diracs);
%Ajouter un vectuer de temps dans l'exponentielle
x = real(xe .* exp(2*1i*pi*Fp*t));
