%% Nettoyage
close all;
clear;

%% Variables initiales
nb_bits = 100;
info_binaire = randi([0,1], 1,nb_bits);
Fe = 24000;
Rb = 3000;
N = 101;

%% Modulateur

% Variables
Ns = Fe/Rb;
a_1 = 1;
a_0 = -1;
h = ones(1,Ns);

% Calculs
mapping = info_binaire.*(a_1 - a_0) + a_0;
Suite_diracs = kron(mapping, [1 zeros(1, Ns-1)]);
x = filter(h, 1, Suite_diracs);
    %DSP
mod_DSP = fftshift(abs(fft(xcorr(x,'unbiased'),1024)));
plage_module=(-Fe/2:Fe/(length(mod_DSP)-1):Fe/2);

% Affichage
figure('Name',"Modulateur");

Bit = [0;1];
Mapping = [a_0;a_1];
T = table(Bit,Mapping);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'FontUnits','normalized','FontSize',0.2,'Units','normalized','Position',[0 0.5 0.5 0.5]);

subplot(2,2,2);
plot((0:1/Fe:(Ns-1)/Fe), h);
title('Filtre de mise en forme rectangulaire');
xlabel('Temps (s)');
ylabel('Hauteur');

subplot(2,2,3)
plot((0:1/Fe:(Ns*nb_bits-1)/Fe), x);
title('Filtrage du modulateur 1');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(2,2,4);
plage_module=(-Fe/2:Fe/(length(mod_DSP)-1):Fe/2);
semilogy(plage_module,mod_DSP);
title("DSP du modulateur");
xlabel('Hz');
ylabel('Module TFD');

%% Demodulateur

%Filtrage de réception
hr = ones(1,Ns);
z = filter(hr, 1, x);

%Echantillonage
z_echant = z(1:Ns:end);

%Décision
info_bin_rec = z_echant > 0;
difference = sum(info_binaire - info_bin_rec);

%Réponse impulsionnelle
g = conv(h,hr);

oeil = reshape(z, 2*Ns, length(z)/(2*Ns));

%Affichage
figure('Name',"Sans canal");
subplot(2,2,[1 2]);
plot((0:1/Fe:(Ns*nb_bits-1)/Fe),z);
title("Signal en sortie de filtre de réception");
xlabel("Temps")
ylabel("Amplitude");
title("Signal filtré");

subplot(2,2,3);
plot((0:1/Fe:2*(Ns-1)/Fe),g);
title('Réponse impulsionnel modulateur+démodulateur');

subplot(2,2,4);
plot(oeil);
title("Diagramme de l'oeil");

n0_1 = Ns;
z_echant_1 = z(n0_1:Ns:end);
info_bin_rec_1 = z_echant_1 > 0;
taux_erreur_binaire_1 = sum(abs(info_bin_rec_1-info_binaire))/length(info_binaire);
fprintf("Taux d'erreur pour n0 = %.1f $%f.\n", n0_1, taux_erreur_binaire_1);

%Test avec n0 = 3

n0_2 = 3;
z_echant_2 = z(n0_2:Ns:end);
info_bin_rec_2 = z_echant_2 > 0;
taux_erreur_binaire_2 = sum(abs(info_bin_rec_2-info_binaire))/length(info_binaire);
fprintf("Taux d'erreur pour n0 = %.1f $%f.\n", n0_2, taux_erreur_binaire_2);

%% Canal de propagation sans bruit

%Calculs communs aux deux sous parties
H = fft(xcorr(h,'unbiased'),1024);
plage_module_HH=(-Fe/2:Fe/(length(H)-1):Fe/2);
Hr = fft(xcorr(hr,'unbiased'),1024);
n0_canal = floor(Ns/2)+floor(Ns/2)+floor(N/2);

% BW = 8000 Hz

fc_1 = 8000;

hc_1 = (2*fc_1/Fe)*sinc(2*(fc_1/Fe)*[-(N-1)/2:(N-1)/2]);
rep_impul_glob_1 = conv(g,hc_1);

    % Passage dans le canal
y_1 = filter(hc_1, 1, [x zeros(1,ceil(n0_canal/(2*Ns))*2*Ns)]); %ajout de suffisament de zéros pour le retard n0 mais aussi suffisament pour être divisible par 2*Ns
z_canal_1 = filter(hr, 1, y_1);
oeil_canal_1 = reshape(z_canal_1, 2*Ns, length(z_canal_1)/(2*Ns));

    % Calcul de la réponse en fréquentielles
Hc_1 = fftshift(abs(fft(xcorr(hc_1,'unbiased'),1024)));
plage_module_H_1=(-Fe/2:Fe/(length(Hc_1)-1):Fe/2);

figure('Name','Passage par canal 8000 Hz');
subplot(2,2,1);
plot(rep_impul_glob_1);
title('Réponse impulsionnelle globale');
subplot(2,2,2);
plot(oeil_canal_1);
title("Diagramme de l'oeil");

subplot(2,2,[3 4]);
sH_1 = semilogy(plage_module_H_1,Hc_1);
hold on;
sHH_1 = semilogy(plage_module_HH,fftshift(abs(Hr.*H)));
hold off;
legend([sH_1, sHH_1],"Réponse du canal","Réponse du modulateur+démodulateur");
title("Réponses fréquentielles");

z_echant_3 = z_canal_1(n0_canal:Ns:end);
info_bin_rec_3 = z_echant_3 > 0;
taux_erreur_binaire_3 = sum(abs(info_bin_rec_3(1:length(info_binaire))-info_binaire))/length(info_binaire);
fprintf("Taux d'erreur pour n0 = %.1f avec un canal filtrant à %.1f Hz : $%f.\n", n0_canal, fc_1, taux_erreur_binaire_3);

% BW = 1000 Hz
fc_2 = 1000;

hc_2 = (2*fc_2/Fe)*sinc(2*(fc_2/Fe)*[-(N-1)/2:(N-1)/2]);
rep_impul_glob_2 = conv(g,hc_2);

    % Passage dans le canal
y_2 = filter(hc_2, 1, [x zeros(1,ceil(n0_canal/(2*Ns))*2*Ns)]);
z_canal_2 = filter(hr, 1, y_2);
oeil_canal_2 = reshape(z_canal_2, 2*Ns, length(z_canal_2)/(2*Ns));

    % Calcul de la réponse en fréquentielles
Hc_2 = fftshift(abs(fft(xcorr(hc_2,'unbiased'),1024)));
plage_module_H_2=(-Fe/2:Fe/(length(Hc_2)-1):Fe/2);

figure('Name','Passage par canal 1000 Hz');
subplot(2,2,1);
plot(rep_impul_glob_2);
title('Réponse impulsionnelle globale');
subplot(2,2,2);
plot(oeil_canal_2);
title("Diagramme de l'oeil");

subplot(2,2,[3 4]);
sH_2 = semilogy(plage_module_H_2,Hc_2);
hold on;
sHH_2 = semilogy(plage_module_HH,fftshift(abs(Hr.*H)));
hold off;
legend([sH_2, sHH_2],"Réponse du canal","Réponse du modulateur+démodulateur");
title("Réponses fréquentielles");

z_echant_4 = z_canal_2(n0_canal:Ns:end);
info_bin_rec_4 = z_echant_4 > 0;
taux_erreur_binaire_4 = sum(abs(info_bin_rec_4(1:length(info_binaire))-info_binaire))/length(info_binaire);
fprintf("Taux d'erreur pour n0 = %.1f avec un canal filtrant à %.1f Hz : $%f.\n", n0_canal, fc_2, taux_erreur_binaire_4);
%pour 1000, les hautes fréquences ne passent pas car filtré par le canal
%(filtrage trop important par rapport au filtrage modulateur+démodulateur)
%donc l'oeil "s'arrondi", on perd de l'information