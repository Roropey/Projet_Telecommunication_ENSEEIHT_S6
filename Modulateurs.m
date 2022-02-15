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
a_1_0 = -10;
a_1_1 = 10;
h_1 = ones(1,Ns_1);

% Calculs
mapping_1 = info_binaire.*(a_1_1 - a_1_0) + a_1_0;
Suite_diracs_1 = kron(mapping_1, [1 zeros(1, Ns_1-1)]);
Suite_diracs_1_decale=[Suite_diracs_1 zeros(1,floor(Ns_1/2))]; 
x_1_decale = filter(h_1, 1, Suite_diracs_1_decale);
%x_1_bis=filter(h_1, 1, Suite_diracs_1);
x_1=x_1_decale(floor(Ns_1/2)+1:end);
%x_1(1:750)-x_1_bis(51:800)
mod_1_DSP = fftshift(abs(fft(xcorr(x_1,'unbiased'))));

% Affichage
figure('Name',"Modulateur 1");

Bit = [0;1];
Mapping = [a_1_0;a_1_1];
T = table(Bit,Mapping);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'FontUnits','normalized','FontSize',0.2,'Units','normalized','Position',[0 0.5 0.5 0.5]);

subplot(2,2,2);
plot(h_1);
title('Filtre de mise en forme rectangulaire');
xlabel('Temps (s)');
ylabel('Hauteur');

subplot(2,2,3)
plot(x_1);
title('Filtrage du modulateur 1');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(2,2,4);
plage_module_1=(-Fe/2:Fe/(length(mod_1_DSP)-1):Fe/2);
semilogy(plage_module_1,mod_1_DSP);
title("DSP du modulateur 1");
xlabel('Hz');
ylabel('Module TFD');

%% Modulateur 2

% Variables
Ns_2 = (Fe/Rb)*2;
a_2_00 =  -10;
a_2_01 = -20;
a_2_10 = 10;
a_2_11 = 20;
h_2 = ones(1,Ns_2);

% Calculs
info_binaire_2 = reshape(info_binaire, [2 nb_bits/2]);
mapping_2 = info_binaire_2(1,:).*info_binaire_2(2,:).*(a_2_11-a_2_10-a_2_01+a_2_00) + info_binaire_2(1,:).*(a_2_10-a_2_00) + info_binaire_2(2,:).*(a_2_01-a_2_00) + a_2_00;
Suite_diracs_2 = kron(mapping_2, [1 zeros(1, Ns_2-1)]);
Suite_diracs_2_decale=[Suite_diracs_2 zeros(1,floor(Ns_2/2))]; 
x_2_decale = filter(h_2, 1, Suite_diracs_2_decale);
x_2=x_2_decale(floor(Ns_2/2)+1:end);
mod_2_DSP = fftshift(abs(fft(xcorr(x_2,'unbiased'))));

% Affichage
figure('Name',"Modulateur 2");

Bit = {'00', '01', '10', '11'}';
Mapping = {int2str(a_2_00),int2str(a_2_01),int2str(a_2_10),int2str(a_2_11)}';
T = table(Bit,Mapping);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'FontUnits','normalized','FontSize',0.15,'Units','normalized','Position',[0 0.5 0.5 0.5]);

subplot(2,2,2);
plot(h_2);
title('Filtre de mise en forme rectangulaire');
xlabel('Temps (s)');
ylabel('Hauteur');

subplot(2,2,3);
plot(x_2);
title('Filtrage du modulateur 2');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(2,2,4);
plage_module_2=(-Fe/2:Fe/(length(mod_2_DSP)-1):Fe/2);
semilogy(plage_module_2,mod_2_DSP);
title("DSP du modulateur 2");
xlabel('Hz');
ylabel('Module TFD');

%% Modulateur 3

% Variables
Ns_3 = Fe/Rb;
a_3_0 = -20;
a_3_1 = 20;
alpha = 0;
h_3 = rcosdesign(alpha, (N-1)/Ns_3,Ns_3);

% Calculs
mapping_3 = info_binaire.*(a_3_1 - a_3_0) + a_3_0;
Suite_diracs_3 = kron(mapping_1, [1 zeros(1, Ns_3-1)]);
Suite_diracs_3_decale=[Suite_diracs_3 zeros(1,floor(N/2))]; 
x_3_decale = filter(h_3, 1, Suite_diracs_3_decale);
x_3=x_3_decale(floor(N/2)+1:end);
mod_3_DSP = fftshift(abs(fft(xcorr(x_3,'unbiased'))));

% Affichage
figure('Name',"Modulateur 3");

Bit = [0;1];
Mapping = [a_3_0;a_3_1];
T = table(Bit,Mapping);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'FontUnits','normalized','FontSize',0.2,'ColumnWidth','auto','Units','normalized','Position',[0 0.5 0.5 0.5]);

subplot(2,2,2);
plot(h_3);
title('Filtre de mise en forme racine de cosinus surélevé');
xlabel('Temps (s)');
ylabel('Hauteur');

subplot(2,2,3);
plot(x_3);
title('Filtrage du modulateur 3');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(2,2,4);
plage_module_3=(-Fe/2:Fe/(length(mod_3_DSP)-1):Fe/2);
semilogy(plage_module_3,mod_2_DSP);
title("DSP du modulateur 3");
xlabel('Hz');
ylabel('Module TFD');

%% Comparaison

% Affichage

figure('Name','Comparaison des modulateurs');
subplot(1,2,1);
mod1 = plot(x_1);
hold on;
mod2 = plot(x_2,'r','Linewidth',1);
mod3 = plot(x_3,'Color',[0.4660 0.6740 0.1880],'Linewidth',1);
hold off;
title("Signaux en sortie des modulateurs")
legend([mod1, mod2, mod3],"Module du modulateur 1","Module du modulateur 2","Module du modulateur 3");
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(1,2,2);
s1 = semilogy(plage_module_1,mod_1_DSP);
hold on;   
s2 = semilogy(plage_module_2,mod_2_DSP,'r','Linewidth',1);
s3 = semilogy(plage_module_3,mod_3_DSP,'Color',[0.4660 0.6740 0.1880],'Linewidth',1);
hold off;
title("DSP des modulateurs");
legend([s1, s2, s3],"Module du modulateur 1","Module du modulateur 2","Module du modulateur 3");
xlabel('Hz');
ylabel('Module TFD');
