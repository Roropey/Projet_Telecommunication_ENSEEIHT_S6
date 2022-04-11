%% Nettoyage
close all;
clear;

%% Récupération information


%try
    fprintf("Execution 5_2\n");
    ChaineDeReference_5_2;
    fprintf("Execution 5_3\n");
    ChaineDeReference_5_3;
    close all;
    %load Chaine_5_2;
    %load Chaine_5_3;
%catch ME
    %if (strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile'))
        %fprintf("Veuillez lancer une première fois les ChaineDeReference_5_x.m .\n");
    %end
%end

%% Affichage de comparaison

figure;
s1 = semilogy(E_bN0dB_2, TEB_5_2);
hold on;
s2 = semilogy(E_bN0dB_3,TEB_5_3);
legend([s1, s2],"Valeur 5_2","Valeur 5_3");
hold off;
xlabel('Eb/N0 (dB)');
ylabel('TEB');
title('TEBs simulés');
