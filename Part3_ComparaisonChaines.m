%% Nettoyage
close all;
clear;

%% Récupération information


try
    load Chaine_5_2;
    load Chaine_5_3;
    load Chaine_5_4_5_6;
catch ME
    if (strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile'))
        printf("Veuillez lancer une première fois les programmes Part3_X.m .\n");
    end
end

%% Affichage de comparaison

figure('Name',"Comparaison de la première étude des taux d'erreur binaires",'Position', [100 100 1300 600]);
s1 = semilogy(E_bN0dB_2, TEB_5_2);
hold on;
s2 = semilogy(E_bN0dB_3,TEB_5_3);
legend([s1, s2],"Valeur de référence","Valeur première étude");
hold off;
xlabel('Eb/N0 (dB)');
ylabel('TEB');
title('TEBs simulés');

figure('Name',"Comparaison de la deuxième étude des taux d'erreur binaires",'Position', [100 100 1300 600]);
s1 = semilogy(E_bN0dB_2, TEB_5_2);
hold on;
s2 = semilogy(E_bN0dB_4,TEB_5_4);
legend([s1, s2],"Valeur de référence","Valeur deuxième étude");
hold off;
xlabel('Eb/N0 (dB)');
ylabel('TEB');
title('TEBs simulés');

figure('Name',"Comparaison des 3 différents taux d'erreur binaires",'Position', [100 100 1300 600]);
s1 = semilogy(E_bN0dB_2, TEB_5_2);
hold on;
s2 = semilogy(E_bN0dB_3,TEB_5_3);
s3 = semilogy(E_bN0dB_4,TEB_5_4);
legend([s1, s2, s3],"Valeur de référence","Valeur première étude","Valeur deuxième étude");
hold off;
xlabel('Eb/N0 (dB)');
ylabel('TEB');
title('TEBs simulés');
