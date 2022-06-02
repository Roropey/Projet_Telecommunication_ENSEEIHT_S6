%% Nettoyage
close all;
clear;

%% Variables initiales
nb_bits = 100;
Fe = 24000;
Rb = 3000;
seuil_derreur = 1000;
a=[-1 1];
Ns = log2(length(a))*Fe/Rb;
h = ones(1,Ns);
hr = h;
n0 = Ns;

%% Calculs
info_binaire = [0 1 1 0 0 1];
[~,info_recu,x,y,z,z_echant] = Propagation_Multi_Canal(info_binaire,nb_bits,Fe,Rb,n0,a,h,hr,1,0,0.5,1/Rb,Inf);

%% Modulation
figure("Name","Sortie modulation");
plot(0:1/Fe:length(info_binaire)/Rb-1/Fe,x);
xlabel("Temps");
ylabel("Amplitude");
title("Signal en sortie de modulation");

%% Canal
figure("Name","Sortie de canal");
plot(0:1/Fe:length(info_binaire)/Rb-1/Fe,y);
xlabel("Temps");
ylabel("Amplitude");
title("Signal en sortie de canal de progation");

%% Demodulateur
figure("Name","Sortie de filtre de récéption");
plot(0:1/Fe:length(info_binaire)/Rb-1/Fe,z);

s=scatterplot(z_echant);
s.Name = "Constellation obtenue en réception";

taux_erreur_binaire = sum(info_recu~=info_binaire)/length(info_binaire);
fprintf("Taux d'erreur pour n0 = %.1f $%f.\n", n0, taux_erreur_binaire);

%% Bruit
TEB = [];
E_bN0_plage = 0:0.1:10;
for E_bN0dB = E_bN0_plage
    nb_bits_faux = 0;
    nb_bits_tot = 0;
    while nb_bits_faux < seuil_derreur
        [info_entree,info_recu,~,~,~,~] = Propagation_Multi_Canal(Inf,nb_bits,Fe,Rb,n0,a,h,hr,1,0,0.5,1/Rb,E_bN0dB);
        nb_bits_faux = nb_bits_faux + sum(info_entree~=info_recu);
        nb_bits_tot = nb_bits_tot + nb_bits;
    end;
    TEB = [TEB nb_bits_faux/nb_bits_tot];
end;

figure;
plot(E_bN0_plage,TEB);



