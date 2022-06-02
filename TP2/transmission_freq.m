function [information_entree information_sortie xe z] = transmission_freq(Fe,Rb,N,type,M,nb_bits,E_bN0Db,n0,h,hr,affichage)
%transmision_freq : réalise la chaîne de transmission d'un message binaire aléatoire sur
%fréquence porteuse
%   Fe : fréquence d'échantillonnage
%   Rb : débit binaire
%   N : ordre de filtrage
%   type : type de modulation en string ('PSK', 'QAM' par exemple)
%   M : ordre de modulation
%   nb_bits : nombre de bits dans le message généré aléatoire
%   E_bN0Db : bruit à appliquer en décibel (Inf pour pas de bruit)
%   n0 : choix échantillonnage
%   h : filtre de mise en forme
%   hr : filtre de réception
%   affichage : booléan pour affichage ou non des constallations (true pour
%   affichage, false pour ne rien afficher)
% -------------------------------------------------------------------------
%   information_entree : message binaire généré aléatoirement
%   information_sortie : message binaire traduit en sortie de chaine de
%   transmission
%   we : signal en sortie de modulation
%   z : signal en sortie de filtre de réception
%   symbole_emis

%% Initialisation

if mod(nb_bits,log2(M))~=0
    error("Nombre de bits incompatible au mappage");
end
information_entree = randi([0,1], 1, nb_bits);
%information_entree_ajout = [zeros(1,log2(M)) information_entree];
%% Modulateur
% Variables
Ns = (Fe/Rb)*log2(M);

info_binaire = reshape(information_entree, [log2(M) (nb_bits)/log2(M)]);

switch type
    case 'ASK'
        if M==4
            mapping = -M/2-1:2:M/2+1;
            symboles_envoye = info_binaire(1,:).*info_binaire(2,:).*(mapping(4)-mapping(3)-mapping(1)+mapping(2))+ info_binaire(1,:).*(mapping(3)-mapping(2)) + info_binaire(2,:).*(mapping(1)-mapping(2)) + mapping(2); 
        else
            error("Que 4-ASK supporter");
        end
    case 'QPSK'
        if M==4
            symboles_envoye = qammod(info_binaire,M,'InputType','bit');
        else
            error("QPSK et M différent de 4 ? problème de définition");
        end
    case 'QAM'
        symboles_envoye = qammod(info_binaire,M,'InputType','bit');
        
    case 'PSK'
        info_dec = bit2int(info_binaire,log2(M));
        symboles_envoye = pskmod(info_dec,M,0,'gray');
        
    otherwise
        error("Type non reconnu ou non supporté : doit être ASK ou QPSK ou QAM ou PSK");
end

if affichage
    info_dec = bit2int(info_binaire,log2(M));
    s=scatterplot(symboles_envoye(1:end));
    % Affichage bit correspondant au mapping /!\ temps de calcul et
    % d'affichage xInf ainsi que forte utilisation de la RAM
%     for k = 1:length(info_dec)-1;
%         switch type
%             case {'QPSK','PSK'}
%                 text(real(symboles_envoye(k))-0.08,imag(symboles_envoye(k))-.08,dec2base(info_dec(k),2,log2(M)),'Color',[1 1 1]);
%             otherwise
%                 text(real(symboles_envoye(k))-0.2,imag(symboles_envoye(k))-.15,dec2base(info_dec(k),2,log2(M)),'Color',[1 1 1]);
%         end
%     end
      s.Name = strcat("Constellation en sortie de mapping : ",int2str(M),'-',type);
end


%% Modulation
Suite_diracs = kron(symboles_envoye, [1 zeros(1, Ns-1)]);
Suite_diracs_decale=[Suite_diracs zeros(1,floor(N/2))]; 
xe_decale = filter(h, 1, Suite_diracs_decale);
xe = xe_decale(floor(N/2)+1:end);

%% Bruit
if E_bN0Db==Inf
    x_bruite = xe;
else
    P_re =  mean(abs(xe).^2);
    Sigma_n = sqrt((P_re*Ns)/(2*log2(M)*10.^(E_bN0Db/10)));
    bruit = Sigma_n*randn(1, length(xe))+1i*Sigma_n*randn(1, length(xe));
    x_bruite = xe + bruit;
end

%% Demodulation
x_bruite_decale = [x_bruite zeros(1,floor(N/2))];
z_decale = filter(hr, 1, x_bruite_decale);
z = z_decale(floor(N/2)+1:end);

z_echant = z(n0:Ns:end);
if affichage
    s=scatterplot(z_echant(2:end));
    s.Name=strcat('Constellation après échantillonnage : ',' ',int2str(M),'-',type,' puissance bruit ',int2str(E_bN0Db));
end
switch type
    case 'ASK'
        symbole_recu = ((abs(real(z_echant))<(mapping(3)+mapping(4))/2)*mapping(3) + (abs(real(z_echant))>=(mapping(3)+mapping(4))/2)*mapping(4)).*sign(real(z_echant));
        information_sortie_reshape=[symbole_recu>(mapping(2)+mapping(3))/2 ; abs(symbole_recu)>(mapping(3)+mapping(4))/2];
    case {'QAM','QPSK'}        
        information_sortie_reshape = qamdemod(z_echant,M,'OutputType','bit');
    case 'PSK'       
       information_sortie_reshape =int2bit(pskdemod(z_echant,M,0,'gray'),log2(M));
end
information_sortie = reshape(information_sortie_reshape,1,nb_bits);
%information_sortie = information_sortie(log2(M)+1:end);
end