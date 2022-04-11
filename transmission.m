function [information_entree information_sortie x y z] = transmission(Fe,Rb,N,a,nb_bits,E_bN0Db,n0,h,h_c,hr)

information_entree = randi([0,1], 1, nb_bits);

% Calculs
mapping = information_entree.*(a(2) - a(1)) + a(1);
Suite_diracs = kron(mapping, [1 zeros(1, Fe/Rb-1)]);
x = filter(h, 1, Suite_diracs);

%% Canal de Transmission
if (E_bN0Db < 0)
    x_bruite = x;
else
    P_x = mean(abs(x).^2);
    Sigma_n = sqrt((P_x*Fe/Rb)/(2*log2(length(a))*10.^(E_bN0Db/10)));
    bruit = Sigma_n*randn(1, length(x));
    x_bruite = x + bruit;
end;
%% Demodulateur

%Filtrage de réception
z = filter(hr, 1, x_bruite);

%Echantillonage

z_echant = z(n0:Fe/Rb:end);
information_sortie = z_echant > 0;   

y=0;
end