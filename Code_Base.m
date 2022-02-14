
nb_bits = 100;

info_binaire = randi([0,1],nb_bits,1);

%binaire
a_1 = -100;
a_0 = 100;

mapping_2 = info_binaire.*(a_1-a_0) + a_0;

%4-aire
a_00 = -10;
a_01 = -100;
a_10 = 10;
a_11 = 100;
info_binaire_4_aire = reshape(info_binaire,[nb_bits/2 2]);
mapping_4 = info_binaire_4_aire(:,1).* info_binaire_4_aire(:,2) .*(a_11 - a_01 - a_10 + a_00) + info_binaire_4_aire(:,1) .*(a_10-a_00) + info_binaire_4_aire(:,2).*(a_01-a_00) + a_00;

