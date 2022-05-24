%% Nettoyage
close all;
clear;

%% Variables initiales
nb_bits = 10000;
seuil_erreur = 1000;
Fe = 24000;
Rb = 3000;
N = 101;
a = [-1 1];
h = ones(1,Fe/Rb);
hr = ones(1,Fe/Rb);
hr(length(hr)/2+1:length(hr)) = 0;
n0 = Fe/Rb;

%% Sans bruit
fprintf("Sans bruit\n");
[info_binaire_env_sans_bruit, info_binaire_recu_sans_bruit, ~, ~, x, ~, z] = transmission(Fe,Rb,N,a,nb_bits,-1,n0,h,0,hr);

mod_DSP = fftshift(abs(fft(xcorr(x,'unbiased'))));
plage_module=(-Fe/2:Fe/(length(mod_DSP)-1):Fe/2);

g = conv(h,hr);

oeil = reshape(z, 2*(Fe/Rb), length(z)/(2*(Fe/Rb)));

taux_erreur_binaire = sum(abs(info_binaire_recu_sans_bruit-info_binaire_env_sans_bruit))/length(info_binaire_env_sans_bruit);

%% Affichage sans bruit

figure('Name',"Modulateur",'Position', [100 100 1300 600]);

Bit = [0;1];
Mapping = a';
T = table(Bit,Mapping);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'FontUnits','normalized','FontSize',0.2,'Units','normalized','Position',[0 0.5 0.5 0.5]);

subplot(2,2,2);
plot((0:((Fe/Rb-1)/Fe)/(length(h)-1):(Fe/Rb-1)/Fe), h);
title('Filtre de mise en forme rectangulaire');
xlabel('Temps (s)');
ylabel('Hauteur');

subplot(2,2,3)
plot((0:1/Fe:((Fe/Rb)*nb_bits-1)/Fe), x);
title('Filtrage du modulateur 1');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(2,2,4);
plage_module=(-Fe/2:Fe/(length(mod_DSP)-1):Fe/2);
semilogy(plage_module,mod_DSP);
title("DSP du modulateur");
xlabel('Hz');
ylabel('Module TFD');

figure('Name','Signal filtré','Position', [100 100 1300 600]);
plot((0:1/Fe:((Fe/Rb)*nb_bits-1)/Fe),z);
title("Signal en sortie de filtre de réception");
xlabel("Temps")
ylabel("Amplitude");

figure('Name','Convolution','Position', [100 100 1300 600]);
plot((0:(2*(Fe/Rb-1)/Fe)/(length(g)-1):2*((Fe/Rb)-1)/Fe),g);
xlabel('Temps');
ylabel('Amplitudes');

figure('Name',"Diagramme de l'oeil",'Position', [100 100 1300 600]);
plot(oeil(:,1:2*(Fe/Rb)));
xlabel('Echantillons');
ylabel('Amplitude');

fprintf("Taux d'erreur pour n0 = %.1f : %.4f.\n", n0, taux_erreur_binaire);

%% Avec bruit
fprintf("Avec bruit\n");

TEB_5_3 = [];
Z_5_3 = [];
E_bN0dB_3 = 0:0.5:8; 
% Calculs
for k=E_bN0dB_3
    nb_bits_faux = 0;
    nb_bits_tot = 0;
    while nb_bits_faux < seuil_erreur
        [info_binaire_env, info_binaire_recu, ~, ~, ~, ~, z] = transmission(Fe,Rb,N,a,nb_bits,k,n0,h,0,hr);
        nb_bits_faux = sum(abs(info_binaire_recu-info_binaire_env)) + nb_bits_faux;
        nb_bits_tot = nb_bits_tot + nb_bits;
    end;
    Z_5_3 = [Z_5_3 z'];
    TEB_5_3 = [TEB_5_3 nb_bits_faux/nb_bits_tot];
end;

%% Théorique

TEB_th = qfunc(sqrt(10.^(E_bN0dB_3/10)));

%% Affichage avec bruit

nb_diagramme_oeil = 6;
E_bN0dB_31_oeil = 4:4:16;
figure('Name',"Différents diagrammes de l'oeil",'Position', [100 100 1300 600]);
i = 1;
for k=E_bN0dB_31_oeil
    [~, ~ , ~, ~, ~, ~, z] = transmission(Fe,Rb,N,a,nb_bits,k,n0,h,0,hr);
    oeil_bruit = reshape(z, 2*(Fe/Rb), length(z)/(2*(Fe/Rb)));
    subplot(2,2,i);
    plot(oeil_bruit(:,1:100*(Fe/Rb)));
    title(strcat("Diagramme de l'oeil pour E_b/N_0", num2str(k), "dB"));
    xlabel("Echantillons");
    ylabel("Amplitude");
    i = i + 1;
end;

%figure('Name',"Différents diagrammes de l'oeil",'Position', [100 100 1300 600]);
%for i=1:nb_diagramme_oeil
    %subplot(floor(sqrt(nb_diagramme_oeil)),ceil(sqrt(nb_diagramme_oeil)),i);
    %oeil_bruit =  reshape(Z_5_3(:,floor(i*size(Z_5_3,2)/nb_diagramme_oeil)), 2*(Fe/Rb), length(Z_5_3(:,i*floor(size(Z_5_3,2)/nb_diagramme_oeil)))/(2*(Fe/Rb)));
    %plot(oeil_bruit(:,1:10*(Fe/Rb))); 
    %title(strcat("Diagramme de l'oeil pour E_b/N_0 ",num2str(E_bN0dB_3(floor(i*size(Z_5_3,2)/nb_diagramme_oeil)))," dB"));
    %xlabel("Echantillons");
    %ylabel("Amplitude");
%end;

figure('Name', "Taux Erreur Binaire",'Position', [100 100 1300 600]);
s1 = semilogy(E_bN0dB_3, TEB_5_3);
hold on;

s2 = semilogy(E_bN0dB_3,TEB_th);

legend([s1, s2],"Valeur pratique","Valeur théorique");
hold off;
xlabel('Eb/N0 (dB)');
ylabel('TEB');
title('TEB simulé et théorique');

clearvars -except TEB_5_3 E_bN0dB_3;
save Chaine_5_3;
