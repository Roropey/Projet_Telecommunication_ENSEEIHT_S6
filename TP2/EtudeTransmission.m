function [TEB,signal_emis] = EtudeTransmission(Fe,Rb,N,type,M,nb_bits,n0,alpha,E_bN0db_TEB,E_bN0db_cons,seuil_erreur,TEB_th)
%EtudeTransmission Summary of this function goes here
%   Fe : fréquence d'échantillonnage
%   Rb : débit binaire
%   N : ordre de filtrage
%   type : type de modulation en string ('PSK', 'QAM' par exemple)
%   M : ordre de modulation
%   nb_bits : nombre de bits dans le message généré aléatoire
%   E_bN0db_TEB : plage de variation de E_b/N0 en dB pour l'étude avec bruit
%   n0 : choix échantillonnage
%   h : filtre de mise en forme
%   hr : filtre de réception
%   TEB_th : TEB théorique
% -------------------------------------------------------------------------
%   TEB : taux d'erreur binaire pour l E_b/N0 variant de 0 à 6 dB
%   signal_emis : signal emis en sortie de modulation sans bruit

%% Initialisation

Ns = (Fe/Rb)*log2(M);
h = rcosdesign(alpha, (N-1)/Ns,Ns);
hr = h;

%% Etude sans bruit
fprintf("Etude sans bruit\n");
[information_entree, information_sortie, signal_emis, z] = transmission_freq(Fe,Rb,N,type,M,nb_bits,Inf,n0,h,hr,true);
taux_erreur_binaire = sum(abs(information_entree-information_sortie))/length(information_entree);


%oeil_angle = reshape(angle(z), 2*log2(M)*(Fe/Rb),[]);

%figure('Name','Oeil angle');
%plot(oeil_angle);

fprintf("Taux d'erreur sans bruit pour n0 = %.1f avec modulation %d-%s : %.4f.\n", n0, M,type,taux_erreur_binaire);
fprintf("Taper sur ENTREE pour continuer (/!\\ cela supprimera les courbes affichées)\n");
pause;
close all;
%% Etude avec bruit

% Affichage constellation
fprintf("Affichage de différentes constallations pour différents bruits\n");
for E_bN0 = E_bN0db_cons
    transmission_freq(Fe,Rb,N,type,M,nb_bits,E_bN0,n0,h,hr,true);
end;

fprintf("Taper sur ENTREE pour continuer (/!\\ cela supprimera les courbes affichées)\n");
pause;
close all;
% Calcul TEB
fprintf("Calculs et affichage du taux d'erreur binaires et symboles\n");
TEB = [];
TES = [];

for E_bN0 = E_bN0db_TEB
    nb_bits_faux = 0;
    nb_bits_tot = 0;
    nb_symboles_faux = 0;
    nb_symboles_tot = 0;
    while nb_bits_faux < seuil_erreur
        [information_entree, information_sortie, ~, ~, symbole_emis, symbole_recu] = transmission_freq(Fe,Rb,N,type,M,nb_bits,E_bN0,n0,h,hr,false);
        
        nb_bits_faux = sum(abs(information_entree-information_sortie)) + nb_bits_faux;
        nb_bits_tot = nb_bits_tot + nb_bits;
        nb_symboles_faux = sum((real(symbole_emis)~=real(symbole_recu))|(imag(symbole_emis)~=imag(symbole_recu))) +nb_symboles_faux;
        nb_symboles_tot = nb_symboles_tot + nb_bits/log2(M);
    end;
    TEB = [TEB nb_bits_faux/nb_bits_tot];
    TES = [TES nb_symboles_faux/nb_symboles_tot];

end;

% Affichage
figure('Name', strcat("Taux Erreur Symbole : ",int2str(M),'-',type));%,'Position', [100 100 1300 600]);
s1_TES = semilogy(E_bN0db_TEB,TES);
hold on;

s2_TES = semilogy(E_bN0db_TEB,TEB_th*log2(M));

hold off;
xlabel('Eb/N0 (dB)');
ylabel('TES');
title('TES simulé et théorique');
legend([s1_TES s2_TES],"Valeur pratique","Valeur théorique");

% Affichage
figure('Name', strcat("Taux Erreur Binaire : ",int2str(M),'-',type));%,'Position', [100 100 1300 600]);
s1_TEB = semilogy(E_bN0db_TEB,TEB);
hold on;

s2_TEB = semilogy(E_bN0db_TEB,TEB_th);

hold off;
xlabel('Eb/N0 (dB)');
ylabel('TEB');
title('TEB simulé et théorique');
legend([s1_TEB s2_TEB],"Valeur pratique","Valeur théorique");

end

