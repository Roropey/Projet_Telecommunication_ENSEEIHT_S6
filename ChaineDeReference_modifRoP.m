%% Nettoyage
close all;
clear;

%% Variables initiales
nb_bits = 10000;
seuil_erreur = 1;
Fe = 24000;
Rb = 3000;
N = 101;
a=[-1 1];

TEB = [];
E_bN0Db = 0:0.1:8;
for k=E_bN0Db
    nb_bits_faux = 0;
    nb_bits_tot = 0;
    while nb_bits_faux < seuil_erreur
        [info_binaire_env info_binaire_recu]=transmission(Fe,Rb,N,a,nb_bits,k);
        nb_bits_faux = sum(abs(info_binaire_env-info_binaire_recu)) + nb_bits_faux;
        nb_bits_tot = nb_bits_tot + nb_bits;
    end;
    fprintf("%f\n",nb_bits_faux);
    TEB = [TEB nb_bits_faux/nb_bits_tot];
end;


figure;
semilogy(E_bN0Db, TEB);



function sigma_n = Sigma_n()
