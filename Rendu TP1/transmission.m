function [information_entree information_sortie symboles_envoye symbole_recu x z] = transmission(Fe,Rb,N,mapping,nb_bits,E_bN0Db,n0,h,h_c,hr)

% Génération d'un message aléatoire
information_entree = randi([0,1], 1, nb_bits);

% Mapping
if (length(mapping)==2)
    symboles_envoye = information_entree.*(mapping(2) - mapping(1)) + mapping(1);
    Ns = Fe/Rb;
elseif (length(mapping)==4)
    information_entree_reshape = reshape(information_entree, [2 nb_bits/2]);
    symboles_envoye = information_entree_reshape(1,:).*information_entree_reshape(2,:).*(mapping(4)-mapping(3)-mapping(1)+mapping(2))+ information_entree_reshape(1,:).*(mapping(3)-mapping(2)) + information_entree_reshape(2,:).*(mapping(1)-mapping(2)) + mapping(2); 
    Ns = 2*Fe/Rb;
end;
% Sur échantillonnage + ajout de zéros
Suite_diracs = kron(symboles_envoye, [1 zeros(1, Ns-1)]);
% Filtre de mise en forme
x = filter(h, 1, Suite_diracs);

%% Canal de Transmission
if (E_bN0Db == Inf)
    % Bas de bruit
    x_bruite = x;
else
    % Calcul puissance signal
    P_x = mean(abs(x).^2);
    % Calcul puissance bruit
    Sigma_n = sqrt((P_x*Ns)/(2*log2(length(mapping))*10.^(E_bN0Db/10)));
    % Calcul bruit
    bruit = Sigma_n*randn(1, length(x));
    % Ajout bruit
    x_bruite = x + bruit;
end;
%% Demodulateur

%Filtrage de réception
z = filter(hr, 1, x_bruite);

%Echantillonage
facteur = n0;
z_echant = z(n0:Ns:end);
% Décision + démapping
if (length(mapping)==2)
    symbole_recu = (z_echant/facteur > (mapping(1)+mapping(2))/2)*mapping(2) + (z_echant/facteur <= (mapping(1)+mapping(2))/2)*mapping(1);
    information_sortie = (symbole_recu == mapping(2));
elseif (length(mapping)==4)
    symbole_recu = ((abs(z_echant/facteur)<(mapping(3)+mapping(4))/2)*mapping(3) + (abs(z_echant/facteur)>=(mapping(3)+mapping(4))/2)*mapping(4)).*sign(z_echant);
    information_sortie_reshape=[symbole_recu>(mapping(2)+mapping(3))/2 ; abs(symbole_recu)>(mapping(3)+mapping(4))/2];
    %information_sortie_reshape = [z_echant/facteur>(mapping(2)+mapping(3))/2 ; abs(z_echant/facteur)>(mapping(3)+mapping(4))/2];
    information_sortie = reshape(information_sortie_reshape,1,nb_bits);
end;
end