%% Nettoyage
close all;
clear;

%% Variables initiales
nb_bits = 100000;
info_binaire = randi([0,1], 1,nb_bits);
Fe = 24000;
Rb = 3000;
N = 101;
Tb = 1/Rb;

%% Modulateur

% Variables
Ns = Fe/Rb;
a_1 = 1;
a_0 = -1;
h = ones(1,Ns);
M = 2;

% Calculs
mapping = info_binaire.*(a_1 - a_0) + a_0;
Suite_diracs = kron(mapping, [1 zeros(1, Ns-1)]);
x = filter(h, 1, Suite_diracs);
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

%% Canal de Transmission

%Bruit
N0 = 0.001;
P_x = mean(abs(x).^2);
E_b = P_x*Tb;




E_bN0Db = 10*log10(E_b./N0);

%Calcul TEB simulé
Eb_Db = 0:0.1:8;
TEB = [];

for i = Eb_Db
    Sigma_n = sqrt((P_x*Ns)/(2*log2(M)*10.^(i/10)));
    bruit = Sigma_n*randn(1, length(x));
    x_bruite = x + bruit;
    hr = ones(1,Ns);
    z = filter(hr, 1, x_bruite);
    n0 = Ns;
    z_echant = z(n0:Ns:end);
    info_bin_rec = z_echant > 0;
    TEB = [TEB sum(abs(info_bin_rec-info_binaire))/length(info_binaire)];
end;
figure;
semilogy(Eb_Db, TEB);
title('TEB simulé');
%% Demodulateur

%Filtrage de réception
hr = ones(1,Ns);
z = filter(hr, 1, x_bruite);

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
fprintf("Taux d'erreur pour n0 = %.1f $%f.\n", n0, taux_erreur_binaire_1);


