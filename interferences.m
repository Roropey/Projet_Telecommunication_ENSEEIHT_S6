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
plot(h);
title('Filtre de mise en forme rectangulaire');
xlabel('Temps (s)');
ylabel('Hauteur');

subplot(2,2,3)
plot(x);
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
plot(z);
title("Signal en sortie de filtre de réception");
xlabel("Temps")
ylabel("Amplitude");

%Echantillonage
plage_echant = [Ns/(4*Fe):Ns/Fe:nb_bits*Ns/Fe];
z_echant = reshape(z,[8, 100]);
z_moy = mean(z_echant,1);

%Décision
info_bin_rec = z_moy > 0;
difference = sum(info_binaire - info_bin_rec);

%Réponse impulsionnelle
g = conv(h,hr);
figure();
plot(g);

oeil = reshape(z,length(z)/Ns, Ns)
figure();
plot(oeil);
