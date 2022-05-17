%% Nettoyage
close all;
clear;

%% Variables Initiales
alpha  = 0.35;
Fp = 2000;
Fe = 10000;
Rb = 2000;
nb_bits = 10000;
%info_binaire = zeros(1,nb_bits);
%info_binaire(10)=1;
info_binaire = randi([0,1], 1,nb_bits);
N = 201;
seuil_erreur = 1000;




[information_entree information_sortie xe z] = transmission_freq(Fe,Rb,N,type,M,nb_bits,E_bN0Db,n0,h,hr,alpha)
