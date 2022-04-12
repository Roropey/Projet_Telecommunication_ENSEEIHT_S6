%% Nettoyage
close all;
clear;

%% Variables Initiales
alpha  = 0.35;
Fp = 2000;
Fe = 10000;
Rb = 2000;
nb_bits = 10000;
info_binaire = randi([0,1], 1,nb_bits);
N = 101;

%% Modulateur

% Variables
Ns = (Fe/Rb)*2;

% 00
a_00 = -1;
b_00 = -1;

% 01
a_01 = -1;
b_01 = 1;

% 11
a_11 = 1;
b_11 = 1;

% 10
a_10 = 1;
b_10 = -1;

h = rcosdesign(alpha, (N-1)/Ns,Ns);
M = 4;
%t = [0:1/Fe:1/Fe*nb_bits/2-1];

% Modulation sur fréquence porteuse
info_binaire_2 = reshape(info_binaire, [2 nb_bits/2]);
mapping = (info_binaire_2(1, :).* (a_11 - a_01) + a_01) + 1i*(info_binaire_2(2, :).* (b_11 - b_10) + b_10);
Suite_diracs = kron(mapping, [1 zeros(1, Ns-1)]);
xe = filter(h, 1, Suite_diracs);

t = (0:length(xe) - 1) / Fe;

I = real(xe);
Q = imag(xe);

x = real(xe .* exp(2*1i*pi*Fp*t));

% DSP
DSP = fftshift(abs(fft(xcorr(x,'unbiased'),10000)));
plage=(-Fe/2 : Fe/2 - 1) * Fe/(length(DSP)-1);

% Retour en bande de base
xcos = x.*cos(2*pi*Fp*t);
xsin = x.*sin(2*pi*Fp*t);
    % Faire un filtrage passe-bas
    % Multiplier xsin par -j et faire la somme

% Démodulation

%syms expr_th_1(f);
%expr_th_1(f) = var(mapping_1)*(Ns_1/Fe).*(sinc(f *(Ns_1/Fe))).^2;

% Affichage
figure('Name', 'Signal modulé')

subplot(2,2,1);
plot(t,I);
title('I(t)');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(2,2,2);
plot(t, Q)
title('Q(t)');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(2,2,3);
plot(t,x);
title('Signal modulé sur porteuse en temporel');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(2,2,4);
semilogy(plage, DSP);
