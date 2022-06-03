%% Nettoyage
close all;
clear;

%% Récupération données
try
    load SurFreqPorteuse;
    load PasseBasEquivalent;
catch ME
    if (strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile'))
        printf("Veuillez lancer une première fois les programmes ChaineSurFreqPort(Eq).m .\n");
    end
end

%% Affichage TEB

figure('Name', "Taux Erreur Binaire sur fréquence et équivalence",'Position', [100 100 1300 600]);
s1_TEB = semilogy(E_bN0db_TEB_eq,TEB_eq);
hold on;

s2_TEB = semilogy(E_bN0db_TEB_sur_freq,TEB_sur_freq);

hold off;
xlabel('Eb/N0 (dB)');
ylabel('TEB');
title('TEB sur fréquence porteuse et sa version équivalente passe-bas');
legend([s1_TEB s2_TEB],"Valeur sur fréquence porteuse","Valeur passe-bas équivalent");