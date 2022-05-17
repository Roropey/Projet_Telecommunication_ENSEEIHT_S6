function [information_entree information_sortie xe z] = transmission_freq(Fe,Rb,N,type,M,nb_bits,E_bN0Db,n0,h,hr,alpha,affichage)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if mod(nb_bits,log2(M))~=0
    error("Nombre de bits incompatible au mappage");
end
information_entree = randi([0,1], 1, nb_bits);
%% Modulateur
% Variables
Ns = (Fe/Rb)*log2(M);

info_binaire = reshape(information_entree, [log2(M) nb_bits/log2(M)]);

if (type == 'A') & (M==4)
    mapping = -M/2+1:2:M/2+1;
    symboles_envoye = information_entree_reshape(1,:).*information_entree_reshape(2,:).*(mapping(4)-mapping(3)-mapping(1)+mapping(2))+ information_entree_reshape(1,:).*(mapping(3)-mapping(2)) + information_entree_reshape(2,:).*(mapping(1)-mapping(2)) + mapping(2); 
    

end


%% Modulation

symboles_envoye = (info_binaire_2(1, :).* (a_11 - a_01) + a_01) + 1i*(info_binaire_2(2, :).* (b_11 - b_10) + b_10);
Suite_diracs = kron(mapping, [1 zeros(1, Ns-1)]);
Suite_diracs_decale=[Suite_diracs zeros(1,floor(N/2))]; 
xe_decale = filter(h, 1, Suite_diracs_decale);
xe = xe_decale(floor(N/2)+1:end);

%% Bruit
if E_bN0Db < 0
    x_bruite = xe;
else
    P_re =  mean(abs(xe).^2);
    Sigma_n = sqrt((P_re*2*Fe/Rb)/(2*log2(M)*10.^(E_bN0/10)));
    bruit = Sigma_n*randn(1, length(x))+1i*Sigma_n*randn(1, length(x));
    x_bruite = xe + bruit;
end

%% Demodulation
x_bruite_decale = [x_bruite zeros(1,floor(N/2))];
z_decale = filter(hr, 1, x_bruite_decale);
z = z_decale(floor(N/2)+1:end);%

z_echant = z(n0:Ns:end);

if (type=='A') & (M==4)
    symbole_recu = ((abs(z_echant/facteur)<(mapping(3)+mapping(4))/2)*mapping(3) + (abs(z_echant/facteur)>=(mapping(3)+mapping(4))/2)*mapping(4)).*sign(z_echant);
    information_sortie_reshape=[symbole_recu>(mapping(2)+mapping(3))/2 ; abs(symbole_recu)>(mapping(3)+mapping(4))/2];
    %information_sortie_reshape = [z_echant/facteur>(mapping(2)+mapping(3))/2 ; abs(z_echant/facteur)>(mapping(3)+mapping(4))/2];
    information_sortie = reshape(information_sortie_reshape,1,nb_bits);
elseif ((type=='Q') | (type=='P')) & (M==4)
    z_reel_recu = real(z_echant) > 0;
    z_imag_recu = imag(z_echant) > 0;
    z_recu = [z_reel_recu; z_imag_recu];
    information_sortie = reshape(z_recu, 1, nb_bits);
end
end