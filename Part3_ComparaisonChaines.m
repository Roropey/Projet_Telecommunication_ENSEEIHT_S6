%% Nettoyage
close all;
clear;

%% Récupération information


%try
    fprintf("Execution de la chaîne de référence\n");
    Part3_ChaineDeReference;
    fprintf("Execution de la première chaîne étudiée\n");
    Part3_PremiereChaineEtudiee;
    fprintf("Execution de la seconde chaîne étudiée\n");
    Part3_DeuxiemeChaineEtudiee;
    close all;
    %load Chaine_5_2;
    %load Chaine_5_3;
%catch ME
    %if (strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile'))
        %fprintf("Veuillez lancer une première fois les ChaineDeReference_5_x.m .\n");
    %end
%end

%% Affichage de comparaison

figure('Name',"Comparaison des différents taux d'erreur binaires",'Position', [100 100 1300 600]);
s1 = semilogy(E_bN0dB_2, TEB_5_2);
hold on;
s2 = semilogy(E_bN0dB_3,TEB_5_3);
s3 = semilogy(E_bN0dB_4,TEB_5_4);
legend([s1, s2, s3],"Valeur de référence","Valeur première étude","Valeur deuxième étude");
hold off;
xlabel('Eb/N0 (dB)');
ylabel('TEB');
title('TEBs simulés');
