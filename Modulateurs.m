%% Nettoyage
close all;
clear;

%% Variables initiales
nb_bits = 100;
info_binaire = randi([0,1], 1,nb_bits);
Fe = 24000;
Rb = 3000;

%% Modulateur 1

% Variables
Ns_1 = 10;
a_1_0 = -10;
a_1_1 = 10;
h = ones(1,Ns_1);

% Calculs
mapping_1 = info_binaire.*(a_1_1 - a_1_0) + a_1_0;
Suite_diracs_1 = kron(mapping_1, [1 zeros(1, Ns_1-1)]);
x = filter(h, 1, Suite_diracs_1);

% Affichage
f1=figure('Name',"Modulateur 1");
Bits = ["0";"1"];
Mapping=[a_1_0;a_1_1];

%T = table(Bits,Mapping);

%uitable('Data', T, 'ColumnName', T.Properties.VariableNames, 'Position', [20 255 200 150]);

%subplot(2,2,2);
%plot(x);
%axis([0 nb_bits*Ns_1 -15 15]);

%% Modulateur 2

% Variables
Ns_2 = 10;
a_2_00 =  -10;
a_2_01 = -100;
a_2_10 = 10;
a_2_11 = 100;
h = ones(1,Ns_2);

% Calculs
info_binaire_2 = reshape(info_binaire, [2 nb_bits/2]);
mapping_2 = info_binaire_2(1,:).*info_binaire_2(2,:).*(a_2_11-a_2_10-a_2_01+a_2_00) + info_binaire_2(1,:).*(a_2_10-a_2_00) + info_binaire_2(2,:).*(a_2_01-a_2_00) + a_2_00;
Suite_diracs_2 = kron(mapping_2, [1 zeros(1, Ns_2-1)]);
x = filter(h, 1, Suite_diracs_2);

% Affichage
figure('Name',"Modulateur 2");
plot(x);
axis([0 nb_bits*Ns_2/2 -150 150]);

%% Modulateur 3

% Variables
Ns_3 = 10;
a_3_0 = -20;
a_3_1 = 20;
h = ones(1,Ns_3);

% Calculs
mapping_3 = info_binaire.*(a_3_1 - a_3_0) + a_3_0;
Suite_diracs_3 = kron(mapping_1, [1 zeros(1, Ns_3-1)]);
x = filter(h, 1, Suite_diracs_3);

% Affichage
figure('Name',"Modulateur 3");
plot(x);
axis([0 nb_bits*Ns_3 -25 25]);
