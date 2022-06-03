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
[information_entree, information_sortie, signal_emis, z, cons_envoye, cons_recu] = transmission_freq(Fe,Rb,N,type,M,nb_bits,Inf,n0,h,hr);
taux_erreur_binaire = sum(abs(information_entree-information_sortie))/length(information_entree);


figure('Name',strcat("Constellation en sortie de mapping : ",int2str(M),'-',type))
scatter(real(cons_envoye),imag(cons_envoye));
xlabel("Partie réel");
ylabel("Partie imaginnaire");
title(strcat("Constellation en sortie de mapping : ",int2str(M),'-',type));

%  s.Name = strcat("Constellation en sortie de mapping : ",int2str(M),'-',type);

figure('Name',strcat('Constellation après échantillonnage : ',' ',int2str(M),'-',type))
scatter(real(cons_recu),imag(cons_recu));
xlabel("Partie réel");
ylabel("Partie imaginnaire");
title(strcat("Constellation après échantillonnage : ",int2str(M),'-',type));

% s=scatterplot(z_echant(2:end));
% s.Name=strcat('Constellation après échantillonnage : ',' ',int2str(M),'-',type,' puissance bruit ',int2str(E_bN0Db));


fprintf("Taux d'erreur sans bruit pour n0 = %.1f avec modulation %d-%s : %.4f.\n", n0, M,type,taux_erreur_binaire);

%% Etude avec bruit

% Affichage constellation
fprintf("Affichage de différentes constallations pour différents bruits\n");
f1=figure('Name','Constellations en sortie de mapping pour différent bruit');
f2=figure('Name','Constellations après échantillonnage pour différent bruit');
for k=1:length(E_bN0db_cons)
    [~, ~, ~, ~, cons_envoye, cons_recu]=transmission_freq(Fe,Rb,N,type,M,nb_bits,E_bN0db_cons(k),n0,h,hr);
    figure(f1);
    subplot(ceil(sqrt(length(E_bN0db_cons))),ceil(sqrt(length(E_bN0db_cons))),k);
    scatter(real(cons_envoye),imag(cons_envoye));
    xlabel("Partie réel");
    ylabel("Partie imaginnaire");
    title(strcat("Constellation en sortie de mapping pour un bruit de ",int2str(E_bN0db_cons(k)),'dB'));
    figure(f2);
    subplot(ceil(sqrt(length(E_bN0db_cons))),ceil(sqrt(length(E_bN0db_cons))),k);
    scatter(real(cons_recu),imag(cons_recu));
    xlabel("Partie réel");
    ylabel("Partie imaginnaire");
    title(strcat("Constellation après échantillonnage pour un bruit de ",int2str(E_bN0db_cons(k)),'dB'));


end;

% Calcul TEB
fprintf("Calculs et affichage du taux d'erreur binaires et symboles\n");
TEB = [];

for E_bN0 = E_bN0db_TEB
    nb_bits_faux = 0;
    nb_bits_tot = 0;
    while nb_bits_faux < seuil_erreur
        [information_entree, information_sortie, ~, ~] = transmission_freq(Fe,Rb,N,type,M,nb_bits,E_bN0,n0,h,hr);
        
        nb_bits_faux = sum(abs(information_entree-information_sortie)) + nb_bits_faux;
        nb_bits_tot = nb_bits_tot + nb_bits;
    end;
    TEB = [TEB nb_bits_faux/nb_bits_tot];

end;


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

