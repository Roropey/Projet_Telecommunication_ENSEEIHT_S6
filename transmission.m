function [information_entree information_sortie] = transmission(Fe,Rb,N,a,nb_bits,E_bN0Db)

information_entree = randi([0,1], 1, nb_bits);
Ns = Fe/Rb;
h = ones(1,Fe/Rb);

% Calculs
mapping = information_entree.*(a(2) - a(1)) + a(1);
Suite_diracs = kron(mapping, [1 zeros(1, Ns-1)]);
x = filter(h, 1, Suite_diracs);

%% Canal de Transmission

P_x = mean(abs(x).^2);
E_b = P_x/Rb;
N0 = E_b/(10^(E_bN0Db/10));
Sigma_n = sqrt((P_x*Ns)/(2*log2(length(a))*10.^(E_bN0Db/10)));
bruit = Sigma_n*randn(1, length(x));

x_bruite = x + bruit;

%% Demodulateur

%Filtrage de réception
hr = ones(1,Ns);
z = filter(hr, 1, x_bruite);

%Echantillonage

n0 = Ns;
z_echant = z(n0:Ns:end);
information_sortie = z_echant > 0;

end