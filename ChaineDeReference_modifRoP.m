%% Nettoyage
close all;
clear;

%% Variables initiales
nb_bits = 1000000;
seuil_erreur = 100;
Fe = 24000;
Rb = 3000;
N = 101;
a=[-1 1];

TEB = [];
E_bN0Db = 0:0.01:8;

% Calculs
for k=E_bN0Db
    nb_bits_faux = 0;
    nb_bits_tot = 0;
    while nb_bits_faux < seuil_erreur
        [info_binaire_env info_binaire_recu]=transmission(Fe,Rb,N,a,nb_bits,k);
        nb_bits_faux = sum(abs(info_binaire_recu-info_binaire_env)) + nb_bits_faux;
        nb_bits_tot = nb_bits_tot + nb_bits;
    end;
    TEB = [TEB nb_bits_faux/nb_bits_tot];
end;

%% Théorique

%Sigma_n_th = sqrt((Fe/Rb)./(2*log2(length(a))*10.^(E_bN0Db/10)));
%TEB_th_tr = qfunc(Sigma_n_th);

TEB_th = 2*((length(a)-1)/(length(a)*log2(length(a)))).*qfunc(sqrt(((6*log2(length(a)))/(length(a)*length(a)-1)).*10.^(E_bN0Db/10)));

%% Affichage
figure;
s1 = semilogy(E_bN0Db, TEB);
hold on;

s2 = semilogy(E_bN0Db,TEB_th);
%s3 = semilogy(E_bN0Db,TEB_th_tr);

legend([s1, s2],"Valeur pratique","Valeur théorique");
hold off;
xlabel('Eb/N0 (dB)');
ylabel('TEB');
title('TEB simulé');

