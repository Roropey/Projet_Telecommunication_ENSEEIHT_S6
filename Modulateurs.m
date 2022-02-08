nb_bits = 100;
info_binaire = randi([0,1], 1,nb_bits);
Fe = 24000;
Rb = 3000;

%Modulateur 1
Ns_1 = 1;
a_1 = 10;
a_0 = -10;

mapping_1 = info_binaire.*(a_1 - a_0) + a_0;

Suite_diracs = kron(mapping, [1 zeros(1, Ns_1-1)]);
h = ones(1,Ns_1);

x = filter(h, 1, Suite_diracs);

figure; plot(x);
axis([0 nb_bits -15 15]);

%Modulateur 2
Ns_2 = 1;
a_00 =  -10;
a_01 = -100;
a_10 = 10;
a_11 = 100;
info_binaire_2 = reshape(info_binaire, [nb_bits/2 2]);
mapping_petit = info_binaire(:, 1).*info_binaire(:,2).*(a_11-a_10-a_01+a_00) + info_binaire(:,1).*(a_10-a_00) + info_binaire(:, 2).*(a_01-a_00) + a_00;
Suite_diracs = kron(mapping, [1 zeros(1, Ns_2-1)]);
h = ones(1,Ns_2);
x = filter(h, 1, Suite_diracs);
figure; plot(x);
axis([0 nb_bits -150 150]);


