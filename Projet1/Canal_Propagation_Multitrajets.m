%% Nettoyage
close all;
clear;

%% Variables initiales
nb_bits = 100;
info_binaire = randi([0,1], 1,nb_bits);
info_binaire = [0 1 1 0 0 1];
Fe = 24000;
Rb = 3000;
N = 101;

%% Modulateur

% Variables
Ns = Fe/Rb;
Ts = Ns/Fe;
a_1 = 1;
a_0 = -1;
h = ones(1,Ns);

% Calculs
mapping = info_binaire.*(a_1 - a_0) + a_0;
Suite_diracs = kron(mapping, [1 zeros(1, Ns-1)]);
x = filter(h, 1, Suite_diracs);
figure("Name","Sortie modulation");
subplot(1,2,1);
plot(x)
subplot(1,2,2);
plot(h);

%% Canal
alpha_0 = 1
tau_0 = 0;
alpha_1 = 0.5;
tau_1 = tau_0 + Ts;

hc = zeros(1,(max(tau_0,tau_1)+Ts)*Fe);
hc(tau_0*Fe+1:(tau_0+Ts)*Fe+1) = alpha_0;
hc(tau_1*Fe+1:(tau_1+Ts)*Fe+1) = alpha_1;

y=filter(hc,1,x);

figure("Name","Sortie de canal");
g_canal=conv(h,hc);
subplot(2,2,1)
plot(hc);
subplot(2,2,2);
plot(g_canal);
subplot(2,2,[3 4]);
plot(y);

%% Demodulateur

%Filtrage de réception
hr = ones(1,Ns);
z = filter(hr, 1, y);
g = conv(g_canal,hr);
figure("Name","Sortie de filtre de récéption");
subplot(2,2,1)
plot(hr);
subplot(2,2,2);
plot(g);
subplot(2,2,3);
oeil = reshape(z, 4*Ns, []);
plot(oeil);
subplot()
subplot(2,2,4);
plot(z);

n0 = Ns;
z_echant = z(n0:Ns:end);
info_bin_rec = z_echant > 0;
taux_erreur_binaire = sum(abs(info_bin_rec-info_binaire))/length(info_binaire);
fprintf("Taux d'erreur pour n0 = %.1f $%f.\n", n0, taux_erreur_binaire);