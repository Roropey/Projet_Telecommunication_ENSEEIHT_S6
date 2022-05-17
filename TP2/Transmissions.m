%% Nettoyage
close all;
%clear;

%% Variables Initiales
alpha  = 0.5;
Fp = 2000;
Fe = 96000;
Rb = 48000;
nb_bits = 120;
%info_binaire = zeros(1,nb_bits);
%info_binaire(10)=1;
info_binaire = randi([0,1], 1,nb_bits);
N = 201;
seuil_erreur = 1000;
n0 = 1;
Ns = (Fe/Rb)*2;
E_bN0Db = -1;
h = rcosdesign(alpha, (N-1)/Ns,Ns);
hr = h;


%[information_entree information_sortie xe z] = transmission_freq(Fe,Rb,N,'ASK',4,nb_bits,E_bN0Db,n0,h,hr,true);


%taux_erreur_binaire = sum(abs(information_entree-information_sortie))/length(information_entree)

%[information_entree information_sortie xe z] = transmission_freq(Fe,Rb,N,'QPSK',4,nb_bits,E_bN0Db,n0,h,hr,true);


%taux_erreur_binaire = sum(abs(information_entree-information_sortie))/length(information_entree)

%[information_entree information_sortie xe z] = transmission_freq(Fe,Rb,N,'PSK',8,nb_bits,E_bN0Db,n0,h,hr,true);


%taux_erreur_binaire = sum(abs(information_entree-information_sortie))/length(information_entree)

[information_entree information_sortie xe z] = transmission_freq(Fe,Rb,N,'QAM',16,nb_bits,E_bN0Db,1,h,hr,true);

taux_erreur_binaire = sum(abs(information_entree-information_sortie))/length(information_entree)


oeil_reel = reshape(real(z), 8*(Fe/Rb), length(z)/(8*(Fe/Rb)));
figure();
plot(oeil_reel);
oeil_imag = reshape(imag(z), 8*(Fe/Rb), length(z)/(8*(Fe/Rb)));
figure();
plot(oeil_imag);