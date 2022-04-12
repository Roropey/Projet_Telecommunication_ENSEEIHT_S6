%% Nettoyage
close all;
%clear;

%% Variables initiales
nb_bits = 10000;
seuil_erreur = 1000;
Fe = 24000;
Rb = 3000;
N = 101;
a = [-1 1];
h = ones(1,Fe/Rb);
%h = ones(1,N);
hr = ones(1,Fe/Rb);
%hr = ones(1,N);
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

figure('Name',"Modulateur");

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

figure('Name','Signal filtré');
plot((0:1/Fe:((Fe/Rb)*nb_bits-1)/Fe),z);
title("Signal en sortie de filtre de réception");
xlabel("Temps")
ylabel("Amplitude");

figure('Name','Convolution');
plot((0:(2*(Fe/Rb-1)/Fe)/(length(g)-1):2*((Fe/Rb)-1)/Fe),g);

figure('Name',"Diagramme de l'oeil");
plot(oeil);

fprintf("Taux d'erreur pour n0 = %.1f : %.4f.\n", n0, taux_erreur_binaire);

%% Avec bruit
fprintf("Avec bruit\n");

TEB_5_2 = [];
E_bN0dB_2 = 0:0.5:8; 

% Calculs
for k=E_bN0dB_2
    nb_bits_faux = 0;
    nb_bits_tot = 0;
    while nb_bits_faux < seuil_erreur
        [info_binaire_env, info_binaire_recu, ~, ~, x, ~, z] = transmission(Fe,Rb,N,a,nb_bits,k,n0,h,0,hr);
        nb_bits_faux = sum(abs(info_binaire_recu-info_binaire_env)) + nb_bits_faux;
        nb_bits_tot = nb_bits_tot + nb_bits;
    end;
    TEB_5_2 = [TEB_5_2 nb_bits_faux/nb_bits_tot];
end;

%% Théorique

%Sigma_n_th = sqrt((Fe/Rb)./(2*log2(length(a))*10.^(E_bN0dB/10)));
%TEB_th_tr = qfunc(Sigma_n_th);

%TEB_th = 2*((length(a)-1)/(length(a)*log2(length(a)))).*qfunc(sqrt(((6*log2(length(a)))/(length(a)*length(a)-1)).*10.^(E_bN0dB/10)));
TEB_th = qfunc(sqrt(2.*10.^(E_bN0dB_2/10)));

%% Affichage sans bruit
figure;
s1 = semilogy(E_bN0dB_2, TEB_5_2);
hold on;

s2 = semilogy(E_bN0dB_2,TEB_th);
%s3 = semilogy(E_bN0dB,TEB_th_tr);
%s3 = semilogy(E_bN0dB,TEB_th_bis);

legend([s1, s2],"Valeur pratique","Valeur théorique");
hold off;
xlabel('Eb/N0 (dB)');
ylabel('TEB');
title('TEB simulé');

%clearvars -except TEB_5_2 E_bN0dB_2;
%save Chaine_5_2;