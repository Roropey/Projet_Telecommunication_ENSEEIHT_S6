%% Nettoyage
close all;
clear;

%% Variables Initiales
alpha  = 0.35;
Fp = 2000;
Fe = 10000;
Rb = 2000;
nb_bits = 100;
info_binaire = randi([0,1], 1,nb_bits);
N = 101;
seuil_erreur = 100;

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
Suite_diracs_decale=[Suite_diracs zeros(1,floor(N/2))]; 
xe_decale = filter(h, 1, Suite_diracs_decale);
xe=xe_decale(floor(N/2)+1:end);

t = (0:length(xe) - 1) / Fe;

I = real(xe);
Q = imag(xe);

x = real(xe .* exp(2*1i*pi*Fp*t));

% DSP
DSP = fftshift(abs(fft(xcorr(x,'unbiased'),10000)));
plage=(-Fe/2 : Fe/2 - 1) * Fe/(length(DSP)-1);

syms expr_th_xe(f);
expr_th_xe(f) = piecewise( abs(f)<=(1-alpha)*Fe/(2*Ns), (var(mapping)*Fe/Ns).*(Ns/Fe),...
(abs(f)>=(1-alpha)*Fe/(2*Ns)) & (abs(f)<=(1+alpha)*Fe/(2*Ns)),(var(mapping)*Fe/Ns).* (Ns/(2*Fe))*(1+cos( (pi * Ns / (Fe * alpha))*(abs(f)- ((1-alpha)*Fe )/ (2*Ns) ))),...
(abs(f)<(1-alpha)*Fe/(2*Ns)) | (abs(f)>(1+alpha)*Fe/(2*Ns)),0);
syms expr_th_x(f);
expr_th_x(f) = 0.25*( ...
piecewise( abs(f-Fp)<=(1-alpha)*Fe/(2*Ns), (var(mapping)*Fe/Ns).*(Ns/Fe),...
(abs(f-Fp)>=(1-alpha)*Fe/(2*Ns)) & (abs(f-Fp)<=(1+alpha)*Fe/(2*Ns)),(var(mapping)*Fe/Ns).* (Ns/(2*Fe))*(1+cos( (pi * Ns / (Fe * alpha))*(abs(f-Fp)- ((1-alpha)*Fe )/ (2*Ns) ))),...
(abs(f-Fp)<(1-alpha)*Fe/(2*Ns)) | (abs(f-Fp)>(1+alpha)*Fe/(2*Ns)),0) + ...
piecewise( abs(-f-Fp)<=(1-alpha)*Fe/(2*Ns), (var(mapping)*Fe/Ns).*(Ns/Fe),...
(abs(-f-Fp)>=(1-alpha)*Fe/(2*Ns)) & (abs(-f-Fp)<=(1+alpha)*Fe/(2*Ns)),(var(mapping)*Fe/Ns).* (Ns/(2*Fe))*(1+cos( (pi * Ns / (Fe * alpha))*(abs(-f-Fp)- ((1-alpha)*Fe )/ (2*Ns) ))),...
(abs(-f-Fp)<(1-alpha)*Fe/(2*Ns)) | (abs(-f-Fp)>(1+alpha)*Fe/(2*Ns)),0)...
);

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

figure('Name',"comparaison DSP");
s1_3 = semilogy(plage,DSP);
hold on
s2_3 = fplot(expr_th_x, [plage(1) plage(length(plage))]);
set(gca,'YScale','log');
hold off;
legend([s1_3, s2_3],"Valeur pratique","Valeur théorique");
title("DSP");
xlabel('Hz');
ylabel('Module TFD');

% Retour en bande de base
xcos = x.*cos(2*pi*Fp*t);
xsin = x.*sin(2*pi*Fp*t);

% Filtrage passe-bas
Fc = 2200; %3500
Plagef = ((-50:50)'*(1/Fe));
Filt_PB = 2*Fc*sinc(2*Fc*Plagef);
%xcos_filt = filter(Filt_PB, 1, xcos);
%xcos_filt = filter(sinc(Fc*Plagef),1,xcos);
xcos_filt = filtrage(xcos',N,Fc,Fe,"bas",true)';

%xsin_filt = filter(Filt_PB, 1, xsin);
%xsin_filt = filter(sinc(Fc*Plagef),1,xsin);
xsin_filt = filtrage(xsin',N,Fc,Fe,"bas",true)';

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


x_demod = xcos_filt - xsin_filt.*1i;

% Démodulation
hr = rcosdesign(alpha, (N-1)/Ns,Ns);
x_demod_decale = [x_demod zeros(1,floor(N/2))];
z_decale = filter(hr, 1, x_demod_decale);
z=z_decale(floor(N/2)+1:end);

n0 = Ns;
z_echant = z(n0:Ns:end);
z_reel_recu = real(z_echant) > 0;
z_imag_recu = imag(z_echant) > 0;
z_recu = [z_reel_recu; z_imag_recu];
z_recu_reshape = reshape(z_recu, 1, nb_bits);

taux_erreur_binaire = sum(abs(info_binaire-z_recu_reshape))/length(info_binaire);


fprintf("Taux d'erreur pour n0 = %.1f : %.4f.\n", n0, taux_erreur_binaire);


figure('Name',"Comparaison partie réelle");p1=plot(t,I);hold on;p2=plot(t,xcos_filt);hold off; legend([p1, p2],"Signal I","Sortie filtre cos");
figure('Name',"Comparaison partie immaginaire");plot(t,Q);hold on;plot(t,-xsin_filt);hold off;legend([p1, p2],"Signal Q","Sortie filtre sin");