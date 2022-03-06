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
Suite_diracs_decale=[Suite_diracs zeros(1,floor((Ns)/2))]; 
x_decale_mod = filter(h, 1, Suite_diracs_decale);
x=x_decale_mod(floor((Ns)/2)+1:end);
    %DSP
mod_DSP = fftshift(abs(fft(xcorr(x,'unbiased'))));
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
x_decale_demod = [x zeros(1,floor(Ns/2))];
z_decale = filter(hr, 1, x_decale_demod);
z = z_decale(floor(Ns/2)+1:end);

figure('Name','Signal filtré');
plot((0:1/Fe:(Ns*nb_bits-1)/Fe),z);
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
figure('Name','Convolution');
plot((0:1/Fe:2*(Ns-1)/Fe),g);

oeil = reshape(z, 2*Ns, length(z)/(2*Ns));
figure('Name','Oeil de chai plu tro ki');
plot(oeil);

n0 = Ns;
z_echant_1 = z(n0:Ns:end);
info_bin_rec_1 = z_echant_1 > 0;
taux_erreur_binaire_1 = sum(abs(info_bin_rec_1-info_binaire))/length(info_binaire);
fprintf("Taux d'erreur avec décalage avant et après filtrage pour n0 = Ns $%f.\n", taux_erreur_binaire_1);


% n0 est normalement lié généré aux retards par les deux filtres (du modulateur et du
% démodulateur), qui valent chacun Ns/2, donc n0 = Ns, cependant le fait
% d'avoir fait des filtrages où l'on décale de Ns/2 à chaque fois pour
% éviter de perdre de l'information du filtrage, ferait pas en sorte
% d'éviter ce retard dans le n0.
% Cette théorie est renforcée par le fait que choisir n0 = Ns donne un taux
% d'erreur proche de 0,5 (le pire possible), alors que prendre n0 = 1 donne
% un taux d'erreur nul...

%% Test sans réalisé le décalage avant et après chaque filtrage

%% Modulateur

% Calculs
mapping = info_binaire.*(a_1 - a_0) + a_0;
Suite_diracs = kron(mapping, [1 zeros(1, Ns-1)]);
%Suite_diracs_decale=[Suite_diracs zeros(1,floor((Ns)/2))]; 
x = filter(h, 1, Suite_diracs);
%x=x_decale_mod(floor((Ns)/2)+1:end);
    %DSP
mod_DSP = fftshift(abs(fft(xcorr(x,'unbiased'))));
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
%x_decale_demod = [x zeros(1,floor(Ns/2))];
z = filter(hr, 1, x);
%z = z_decale(floor(Ns/2)+1:end);

figure('Name','Signal filtré');
plot((0:1/Fe:(Ns*nb_bits-1)/Fe),z);
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
figure('Name','Convolution');
plot((0:1/Fe:2*(Ns-1)/Fe),g);

oeil = reshape(z, 2*Ns, length(z)/(2*Ns));
figure('Name','Oeil de chai plu tro ki');
plot(oeil);

n0 = Ns;
z_echant_1 = z(n0:Ns:end);
info_bin_rec_1 = z_echant_1 > 0;
taux_erreur_binaire_1 = sum(abs(info_bin_rec_1-info_binaire))/length(info_binaire);
fprintf("Taux d'erreur sans décalage avant et après filtrage pour n0 = Ns $%f.\n", taux_erreur_binaire_1);