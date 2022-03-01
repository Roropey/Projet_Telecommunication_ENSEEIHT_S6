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
mapping_1 = info_binaire.*(a_1 - a_0) + a_0;
Suite_diracs = kron(mapping_1, [1 zeros(1, Ns-1)]);
Suite_diracs_decale=[Suite_diracs zeros(1,floor(Ns/2))]; 
x_decale = filter(h, 1, Suite_diracs_decale);
%x_bis=filter(h, 1, Suite_diracs);
x=x_decale(floor(Ns/2)+1:end);
%x(1:750)-x_bis(51:800)
mod_DSP = fftshift(abs(fft(xcorr(x,'unbiased'))));

% Affichage
figure('Name',"Modulateur");

Bit = [0;1];
Mapping = [a_0;a_1];
T = table(Bit,Mapping);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'FontUnits','normalized','FontSize',0.2,'Units','normalized','Position',[0 0.5 0.5 0.5]);

subplot(2,2,2);
plot((1/Fe:1/Fe:Ns/Fe), h);
title('Filtre de mise en forme rectangulaire');
xlabel('Temps (s)');
ylabel('Hauteur');

subplot(2,2,3)
plot((1/Fe:1/Fe:Ns*nb_bits/Fe), x);
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
x_decale_2 = [x zeros(1,floor(Ns/2))];
z_decale = filter(hr, 1, x_decale_2);
z = z_decale(floor(Ns/2)+1:end);

figure();
plot((1/Fe:1/Fe:Ns*nb_bits/Fe),z);
title("Signal en sortie de filtre de réception");
xlabel("Temps")
ylabel("Amplitude");

%Echantillonage
z_echant = z(1:Ns:end);

%Décision
info_bin_rec = z_echant > 0;
difference = sum(info_binaire - info_bin_rec);

%Réponse impulsionnelle
g = conv(h,hr);
figure();
plot(g);

n0 = Ns;

oeil = reshape(z, 2*Ns, length(z)/(2*Ns));
figure();
plot(oeil);

z_echant_1 = z(Ns:Ns:end);
info_bin_rec_1 = z_echant_1 > 0;
taux_erreur_binaire_1 = sum(abs(info_bin_rec_1-info_binaire))/length(info_binaire)





