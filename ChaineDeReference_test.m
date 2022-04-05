%% Nettoyage
%close all;
%clear;

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
P_x = mean(abs(x).^2)

%Calcul TEB simulé
Eb_Db = 0:0.01:6;
Sigma_n = sqrt((P_x*Ns)./(2*log2(M).*10.^(Eb_Db./10)));
bruit = Sigma_n'*randn(1, length(x));
bruite = repmat(x,size(bruit,1),1) + bruit;
hr = ones(1,Ns);
z = filter(hr, 1, bruite);
n0 = Ns;
z_echant = z(:,n0:Ns:end);
info_bin_rec = z_echant > 0;
TEB_lol = sum(abs(info_bin_rec-repmat(info_binaire,size(info_bin_rec,1),1)),2)./length(info_binaire);
save test;
%figure;
%semilogy(Eb_Db, TEB_lol);
%title('TEB simulé');
