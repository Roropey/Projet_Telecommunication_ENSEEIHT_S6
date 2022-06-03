%% Nettoyage
close all;
clear;

%% Variables initiales
nb_bits = 100;
info_binaire = randi([0,1], 1,nb_bits);
Fe = 24000;
Rb = 3000;
N = 101;

%% Modulateur 1

% Variables
Ns_1 = Fe/Rb;
a_1_0 = -1;
a_1_1 = 1;
h_1 = ones(1,Ns_1);

% Mapping
mapping_1 = info_binaire.*(a_1_1 - a_1_0) + a_1_0;
% Suréchantillonnage
Suite_diracs_1 = kron(mapping_1, [1 zeros(1, Ns_1-1)]);
% Filtre de mise en forme
Suite_diracs_1_decale=[Suite_diracs_1 zeros(1,floor((Ns_1)/2))]; %Ajout de zéros pour ne pas perdre d'information à cause du retard
x_1_decale = filter(h_1, 1, Suite_diracs_1_decale);
x_1=x_1_decale(floor((Ns_1)/2)+1:end);

% DSP pratique
mod_1_DSP = fftshift(abs(fft(xcorr(x_1,'unbiased'),1024)));
plage_module_1=(-Fe/2:Fe/(length(mod_1_DSP)-1):Fe/2);
% DPS théorique
syms expr_th_1(f);
expr_th_1(f) = var(mapping_1)*(Ns_1/Fe).*(sinc(f *(Ns_1/Fe))).^2; 

% Affichage
figure('Name',"Modulateur 1 filtre",'Position', [100 100 1000 600]);

plot((0:1/Fe:(Ns_1-1)/Fe),h_1);
title('Filtre de mise en forme rectangulaire');
xlabel('Temps (s)');
ylabel('Hauteur');

figure('Name',"Modulateur 1 signal",'Position', [100 100 1000 600]);
plot((0:1/Fe:(Ns_1*nb_bits-1)/Fe),x_1);
title('Signal transmis');
xlabel('Temps (s)');
ylabel('Amplitude');

figure('Name',"Modulateur 1 DSP",'Position', [100 100 1000 600]);
semilogy(plage_module_1,mod_1_DSP);
title("DSP du modulateur 1");
xlabel('Hz');
ylabel('Module TFD');

figure('Name',"Modulateur 1 DSP pratique et théorique",'Position', [100 100 1000 600]);
s1_1 = semilogy(plage_module_1,mod_1_DSP);
hold on
s2_1 = fplot(expr_th_1, [plage_module_1(1) plage_module_1(length(plage_module_1))]);
set(gca,'YScale','log');
hold off;
legend([s1_1, s2_1],"Valeur pratique","Valeur théorique");
title("DSP pratique et théorique du modulateur 1");
xlabel('Hz');
ylabel('Module TFD');

%% Modulateur 2

% Variables
Ns_2 = (Fe/Rb)*2;
a_2_00 =  -1;
a_2_01 = -3;
a_2_10 = 1;
a_2_11 = 3;
h_2 = ones(1,Ns_2);


% Mise en forme de l'information binaire en 2 bits par symbole
info_binaire_2 = reshape(info_binaire, [2 nb_bits/2]);
% Mapping
mapping_2 = info_binaire_2(1,:).*info_binaire_2(2,:).*(a_2_11-a_2_10-a_2_01+a_2_00) + info_binaire_2(1,:).*(a_2_10-a_2_00) + info_binaire_2(2,:).*(a_2_01-a_2_00) + a_2_00;
% Sur échantillonnage
Suite_diracs_2 = kron(mapping_2, [1 zeros(1, Ns_2-1)]);
% Filtre de mise en forme
Suite_diracs_2_decale=[Suite_diracs_2 zeros(1,floor((Ns_2)/2))]; %Ajout de zéros pour ne pas perdre d'information à cause du retard
x_2_decale = filter(h_2, 1, Suite_diracs_2_decale);
x_2=x_2_decale(floor((Ns_2)/2)+1:end);
% DSP pratique
mod_2_DSP = fftshift(abs(fft(xcorr(x_2,'unbiased'),1024)));
plage_module_2=(-Fe/2:Fe/(length(mod_2_DSP)-1):Fe/2);
% DSP théorique
syms expr_th_2(f);
expr_th_2(f) = var(mapping_2)*(Ns_2/Fe).*(sinc(f *(Ns_2/Fe))).^2;


% Affichage
figure('Name',"Modulateur 2 filtre",'Position', [100 100 1000 600]);
plot((0:1/Fe:(Ns_2-1)/Fe),h_2);
title('Filtre de mise en forme rectangulaire');
xlabel('Temps (s)');
ylabel('Hauteur');

figure('Name',"Modulateur 2 signal",'Position', [100 100 1000 600]);
plot((0:1/Fe:(Ns_2*nb_bits-1)/(2*Fe)),x_2);
title('Signal transmis');
xlabel('Temps (s)');
ylabel('Amplitude');

figure('Name',"Modulateur 2 DSP",'Position', [100 100 1000 600]);
semilogy(plage_module_2,mod_2_DSP);
title("DSP du modulateur 2");
xlabel('Hz');
ylabel('Module TFD');

figure('Name',"Modulateur 2 DSP pratique et théorique",'Position', [100 100 1000 600]);
s1_2 = semilogy(plage_module_2,mod_2_DSP);
hold on
s2_2 = fplot(expr_th_2, [plage_module_2(1) plage_module_2(length(plage_module_2))]);
set(gca,'YScale','log');
hold off;
legend([s1_2, s2_2],"Valeur pratique","Valeur théorique");
title("DSP théorique et pratique du modulateur 2");
xlabel('Hz');
ylabel('Module TFD');

%% Modulateur 3

% Variables
Ns_3 = Fe/Rb;
a_3_0 = -1;
a_3_1 = 1;
alpha = 0.35;
h_3 = rcosdesign(alpha, (N-1)/Ns_3,Ns_3);


% Mapping
mapping_3 = info_binaire.*(a_3_1 - a_3_0) + a_3_0; 
% Sur échantillonnage
Suite_diracs_3 = kron(mapping_1, [1 zeros(1, Ns_3-1)]);
% Filtre de mise en forme
Suite_diracs_3_decale=[Suite_diracs_3 zeros(1,floor(N/2))]; %Ajout de zéros pour ne pas perdre d'information à cause du retard
x_3_decale = filter(h_3, 1, Suite_diracs_3_decale);
x_3=x_3_decale(floor(N/2)+1:end);
% DSP pratique
mod_3_DSP = fftshift(abs(fft(xcorr(x_3,'unbiased'),1024)));
plage_module_3=(-Fe/2:Fe/(length(mod_3_DSP)-1):Fe/2);
% DSP théorique
syms expr_th_3(f);
expr_th_3(f) = piecewise( abs(f)<=(1-alpha)*Fe/(2*Ns_3), (var(mapping_3)*Fe/Ns_3).*(Ns_3/Fe),...
(abs(f)>=(1-alpha)*Fe/(2*Ns_3)) & (abs(f)<=(1+alpha)*Fe/(2*Ns_3)),(var(mapping_3)*Fe/Ns_3).* (Ns_3/(2*Fe))*(1+cos( (pi * Ns_3 / (Fe * alpha))*(abs(f)- ((1-alpha)*Fe )/ (2*Ns_3) ))),...
(abs(f)<(1-alpha)*Fe/(2*Ns_3)) | (abs(f)>(1+alpha)*Fe/(2*Ns_3)),0);

% Affichage
figure('Name',"Modulateur 3 filtre",'Position', [100 100 1000 600]);

plot((-Ns_3/Fe:2*(Ns_3)/((N-1)*Fe):Ns_3/Fe),h_3);
title('Filtre de mise en forme racine de cosinus surélevé');
xlabel('Temps (s)');
ylabel('Hauteur');

figure('Name',"Modulateur 3 signal",'Position', [100 100 1000 600]);
plot((0:1/Fe:(Ns_3*nb_bits-1)/Fe),x_3);
title('Signal transmis');
xlabel('Temps (s)');
ylabel('Amplitude');

figure('Name',"Modulateur 3 DSP",'Position', [100 100 1000 600]);
semilogy(plage_module_3,mod_3_DSP);
title("DSP du modulateur 3");
xlabel('Hz');
ylabel('Module TFD');

figure('Name',"Modulateur 3 DSP pratique et théorique",'Position', [100 100 1000 600]);
s1_3 = semilogy(plage_module_3,mod_3_DSP);
hold on
s2_3 = fplot(expr_th_3, [plage_module_3(1) plage_module_3(length(plage_module_3))]);
set(gca,'YScale','log');
hold off;
legend([s1_3, s2_3],"Valeur pratique","Valeur théorique");
title("DSP du modulateur 3");
xlabel('Hz');
ylabel('Module TFD');
%% Comparaison

% Affichage

figure('Name','Comparaison des modulateurs signaux','Position', [100 100 1000 600]);
mod1 = plot(x_1);
hold on;
mod2 = plot(x_2,'r','Linewidth',1);
mod3 = plot(x_3,'Color',[0.4660 0.6740 0.1880],'Linewidth',1);
hold off;
title("Signaux en sortie des modulateurs")
legend([mod1, mod2, mod3],"Signal Modulateur 1","Signal Modulateur 2","Signal Modulateur 3");
xlabel('Temps (s)');
ylabel('Amplitude');

figure('Name','Comparaison des modulateurs DSPs','Position', [100 100 1000 600]);
s1 = semilogy(plage_module_1,mod_1_DSP);
hold on;   
s2 = semilogy(plage_module_2,mod_2_DSP,'r','Linewidth',1);
s3 = semilogy(plage_module_3,mod_3_DSP,'Color',[0.4660 0.6740 0.1880],'Linewidth',1);
hold off;
title("DSP des modulateurs");
legend([s1, s2, s3],"DSP Modulateur 1","DSP Modulateur 2","DSP Modulateur 3");
xlabel('Hz');
ylabel('Module TFD');
