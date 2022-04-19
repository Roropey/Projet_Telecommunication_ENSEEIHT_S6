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

% Retour en bande de base
xcos = x.*cos(2*pi*Fp*t);
xsin = x.*sin(2*pi*Fp*t);

% Filtrage passe-bas
Fc = 3500;
Plagef = ((-50:50)'*(1/Fe));
Filt_PB = 2*Fc*sinc(pi*2*Fc*Plagef);
%xcox_filt = filter(Filt_PB, 1, xcos);
xcos_filt = filter(sinc(Fc*Plagef),1,xcos);
%xsin_filt = filter(Filt_PB, 1, xsin);
xsin_filt = filter(sinc(Fc*Plagef),1,xsin);

% Affichage
figure('Name', "Filtrage Passe-Bas");
subplot(2,2,1);
plot(xcos_filt);
title("Xcos filtré, en temporel");
xlabel("Temps");
ylabel("Amplitude");
subplot(2,2,3);
plot(xsin_filt);
title("Xsin filtré, en temporel");
xlabel("Temps");
ylabel("Amplitude");

DSP_cos = fftshift(abs(fft(xcorr(xcos_filt,'unbiased'),10000)));
subplot(2,2,2);
semilogy(plage, DSP_cos);
title("Xcos filtré, en fréquentiel");
DSP_sin = fftshift(abs(fft(xcorr(xsin_filt,'unbiased'),10000)));
subplot(2,2,4);
semilogy(plage, DSP_sin);
title("Xsin filtré, en fréquentiel");

xsin_filt = xsin_filt*1i;

x_demod = xcos_filt - xsin_filt;

% Démodulation
hr = rcosdesign(alpha, (N-1)/Ns,Ns);
z = filter(hr, 1, x_demod);

n0 = 8;
z_echant = z(n0:Ns:end);
z_reel_recu = real(z_echant) > 0;
z_imag_recu = imag(z_echant) > 0;
z_recu = [z_reel_recu; z_imag_recu];
z_recu_reshape = reshape(z_recu, 1, nb_bits);

TEB = [];
E_bN0dB = 0:0.5:8;

for k=E_bN0dB_2
    nb_bits_faux = 0;
    nb_bits_tot = 0;
    while nb_bits_faux < seuil_erreur
        % Mettre les calculs du dessus
        [info_binaire_env, info_binaire_recu, ~, ~, x, ~, z] = transmission(Fe,Rb,N,a,nb_bits,k,n0,h,0,hr);
        nb_bits_faux = sum(abs(info_binaire_recu-info_binaire_env)) + nb_bits_faux;
        nb_bits_tot = nb_bits_tot + nb_bits;
    end;
    TEB = [TEB nb_bits_faux/nb_bits_tot];
end;

%syms expr_th_1(f);
%expr_th_1(f) = var(mapping_1)*(Ns_1/Fe).*(sinc(f *(Ns_1/Fe))).^2;



