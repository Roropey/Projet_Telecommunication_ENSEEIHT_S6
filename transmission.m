function [information_entree information_sortie x y z] = transmission(Fe,Rb,N,a,nb_bits,E_bN0Db,n0,h,h_c,hr)

information_entree = randi([0,1], 1, nb_bits);

% Calculs
if (length(a)==2)
    mapping = information_entree.*(a(2) - a(1)) + a(1);
    Ns = Fe/Rb;
elseif (length(a)==4)
    information_entree_reshape = reshape(information_entree, [2 nb_bits/2]);
    mapping = information_entree_reshape(1,:).*information_entree_reshape(2,:).*(a(4)-a(3)-a(1)+a(2))+ information_entree_reshape(1,:).*(a(3)-a(2)) + information_entree_reshape(2,:).*(a(1)-a(2)) + a(2); 
    Ns = 2*Fe/Rb;
end;
Suite_diracs = kron(mapping, [1 zeros(1, Ns-1)])
x = filter(h, 1, Suite_diracs)

%% Canal de Transmission
if (E_bN0Db < 0)
    x_bruite = x
else
    P_x = mean(abs(x).^2);
    Sigma_n = sqrt((P_x*Ns)/(2*log2(length(a))*10.^(E_bN0Db/10)));
    bruit = Sigma_n*randn(1, length(x));
    x_bruite = x + bruit;
end;
%% Demodulateur

%Filtrage de réception
z = filter(hr, 1, x_bruite);

%Echantillonage

z_echant = z(n0:Ns:end);
if (length(a)==2)
    information_sortie = z_echant > (a(1)+a(2))/2;
elseif (length(a)==4)
    information_sortie_reshape = [z_echant>(a(2)+a(3))/2;abs(z_echant)>8*(a(3)+a(4))/2];
    information_sortie = reshape(information_sortie_reshape,1,nb_bits);
end;
y=0;
end