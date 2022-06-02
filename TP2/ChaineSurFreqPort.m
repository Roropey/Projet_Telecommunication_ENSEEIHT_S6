%% Nettoyage
close all;
clear;

%% Variables Initiales
alpha  = 0.35;
Fp = 2000;
Fe = 10000;
Rb = 2000;
nb_bits = 10000;
%info_binaire = zeros(1,nb_bits);
%info_binaire(10)=1;
info_binaire = randi([0,1], 1,nb_bits);
N = 201;
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

%% Sans bruit

% Modulation sur fréquence porteuse
info_binaire_2 = reshape(info_binaire, [2 nb_bits/2]);
mapping = (info_binaire_2(1, :).* (a_11 - a_01) + a_01) + 1i*(info_binaire_2(2, :).* (b_11 - b_10) + b_10);
Suite_diracs = kron(mapping, [1 zeros(1, Ns-1)]);
Suite_diracs_decale=[Suite_diracs zeros(1,floor(N/2))]; 
xe_decale = filter(h, 1, Suite_diracs_decale);
xe = xe_decale(floor(N/2)+1:end);

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
figure('Name', 'Signal modulé', 'Position', [100 100 1300 600])

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

subplot(2,2,[3,4]);
plot(t,x);
title('Signal modulé sur porteuse en temporel');
xlabel('Temps (s)');
ylabel('Amplitude');

figure('Name', 'DSP du signal modulé sur fréquence porteuse', 'Position', [100 100 1300 600])
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
Filt_PB = 2*Fc*sinc(2*Fc*Plagef)*(1/Fe);
%xcos_filt = filter(Filt_PB, 1, xcos);
xcos_filt = filtrage(xcos',N,Fc,Fe,"bas",true)';

%xsin_filt = filter(Filt_PB, 1, xsin);
xsin_filt = filtrage(xsin',N,Fc,Fe,"bas",false)';

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
hr = h; %rcosdesign(alpha, (N-1)/Ns,Ns);
x_demod_decale = [x_demod zeros(1,floor(N/2))];
z_decale = filter(hr, 1, x_demod_decale);
z = z_decale(floor(N/2)+1:end);%

%oeil_reel = reshape(real(z), 4*(Fe/Rb), length(z)/(4*(Fe/Rb)));
%figure('Name', "Diagramme de l'oeil de la partie réelle du signal en sortie du filtre de réception", 'Position', [100 100 1300 600]);
%plot(oeil_reel(:,1:4*(Fe/Rb)));
%xlabel('Echantillons');
%ylabel('Amplitude');

%oeil_imag = reshape(imag(z), 4*(Fe/Rb), length(z)/(4*(Fe/Rb)));
%figure('Name', "Diagramme de l'oeil de la partie imaginaire du signal en sortie du filtre de réception", 'Position', [100 100 1300 600]);
%plot(oeil_imag(:,1:4*(Fe/Rb)));
%xlabel('Echantillons');
%ylabel('Amplitude');

n0 = 1;
z_echant = z(n0:Ns:end);
z_reel_recu = real(z_echant) > 0;
z_imag_recu = imag(z_echant) > 0;
z_recu = [z_reel_recu; z_imag_recu];
z_recu_reshape = reshape(z_recu, 1, nb_bits);

taux_erreur_binaire = sum(abs(info_binaire-z_recu_reshape))/length(info_binaire);


fprintf("Taux d'erreur sans bruit pour n0 = %.1f : %.4f.\n", n0, taux_erreur_binaire);


figure('Name',"Comparaison partie réelle");p1=plot(t,I);hold on;p2=plot(t,xcos_filt);hold off; legend([p1, p2],"Signal I","Sortie filtre cos");
figure('Name',"Comparaison partie immaginaire");p1=plot(t,Q);hold on;p2=plot(t,-xsin_filt);hold off;legend([p1, p2],"Signal Q","Sortie filtre sin");

%figure();plot(imag(mapping));hold on;plot(imag(z_echant));hold off
%figure();plot(imag(Suite_diracs));hold on;plot(imag(z));hold off

%% Avec bruit



TEB = [];
E_bN0db = 0:0.1:6;


for E_bN0 = E_bN0db
    nb_bits_faux = 0;
    nb_bits_tot = 0;
    while nb_bits_faux < seuil_erreur
        info_binaire = randi([0,1], 1,nb_bits);
        % Modulation sur fréquence porteuse
        info_binaire_2 = reshape(info_binaire, [2 nb_bits/2]);
        mapping = (info_binaire_2(1, :).* (a_11 - a_01) + a_01) + 1i*(info_binaire_2(2, :).* (b_11 - b_10) + b_10);
        Suite_diracs = kron(mapping, [1 zeros(1, Ns-1)]);
        Suite_diracs_decale=[Suite_diracs zeros(1,floor(N/2))]; 
        xe_decale = filter(h, 1, Suite_diracs_decale);
        xe = xe_decale(floor(N/2)+1:end);
        
        x = real(xe .* exp(2*1i*pi*Fp*t));

        P_x =  mean(abs(x).^2);
        Sigma_n = sqrt((P_x*2*Fe/Rb)/(2*log2(M)*10.^(E_bN0/10)));
        bruit = Sigma_n*randn(1, length(x));
        x_bruite = x + bruit;
        %figure();subplot(2,2,1);plot(x);subplot(2,2,2);plot(bruit);subplot(2,2,[3 4]);p1=plot(x);hold on;p2=plot(bruit);hold off;legend([p1,p2],"x","bruit");

        
        % Retour en bande de base
        xcos = x_bruite.*cos(2*pi*Fp*t);
        xsin = x_bruite.*sin(2*pi*Fp*t);

        
        % Filtrage passe-bas
        Fc = 2200;
        xcos_filt = filtrage(xcos',N,Fc,Fe,"bas",false)';
        xsin_filt = filtrage(xsin',N,Fc,Fe,"bas",false)';
        
        
        x_demod = xcos_filt - xsin_filt.*1i;
        
        % Démodulation
        hr = h; %rcosdesign(alpha, (N-1)/Ns,Ns);
        x_demod_decale = [x_demod zeros(1,floor(N/2))];
        z_decale = filter(hr, 1, x_demod_decale);
        z = z_decale(floor(N/2)+1:end);%
        
        n0 = 1;
        z_echant = z(n0:Ns:end);
        z_reel_recu = real(z_echant) > 0;
        z_imag_recu = imag(z_echant) > 0;
        z_recu = [z_reel_recu; z_imag_recu];
        z_recu_reshape = reshape(z_recu, 1, nb_bits);
        
        
        nb_bits_faux = sum(abs(info_binaire-z_recu_reshape)) + nb_bits_faux;
        nb_bits_tot = nb_bits_tot + nb_bits;
    end;
    TEB = [TEB nb_bits_faux/nb_bits_tot];

end;


TEB_th = (4/ log2(M)).*(1-(1/sqrt(M))).*qfunc(sqrt(((3*log2(M))/(M-1)).*10.^(E_bN0db/10)));

figure('Name', "Taux Erreur Binaire",'Position', [100 100 1300 600]);
s1_TEB = semilogy(E_bN0db,TEB);
hold on;

s2_TEB = semilogy(E_bN0db,TEB_th);

hold off;
xlabel('Eb/N0 (dB)');
ylabel('TEB');
title('TEB simulé et théorique');
legend([s1_TEB s2_TEB],"Valeur pratique","Valeur théorique");

