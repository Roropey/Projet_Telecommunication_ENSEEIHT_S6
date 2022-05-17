function [information_entree information_sortie xe z] = transmission_freq(Fe,Rb,N,type,M,nb_bits,E_bN0Db,n0,h,hr,affichage)
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

switch type
    case 'ASK'
        if M==4
            mapping = -M/2-1:2:M/2+1;
            symboles_envoye = info_binaire(1,:).*info_binaire(2,:).*(mapping(4)-mapping(3)-mapping(1)+mapping(2))+ info_binaire(1,:).*(mapping(3)-mapping(2)) + info_binaire(2,:).*(mapping(1)-mapping(2)) + mapping(2); 
            if affichage
                %figure('Name','Constellations en sortie de mapping 4-ASK');
                %plot(symboles_envoye,zeros(1,length(symboles_envoye)),'LineStyle','none','Marker','.');
                scatterplot(symboles_envoye);
                text(mapping(1)-0.12,-0.15,'00','Color','w');
                text(mapping(2)-0.12,-0.15,'01','Color','w');
                text(mapping(3)-0.12,-0.15,'11','Color','w');
                text(mapping(4)-0.12,-0.15,'10','Color','w');
                xlabel("In-phase Amplitude");
                ylabel("Quadrature Amplitude");
            end
        else
            error("Que 4-ASK supporter");
        end
    case 'QPSK'
        if M==4
            symboles_envoye = qammod(info_binaire,M,'InputType','bit');
            %if affichage

              %  for k = 1:M
             %       text(real(symgray(k))-0.2,imag(symgray(k))+.15,...
                     %   dec2base(mapgray(k),2,4));
                    % text(real(symgray(k))-0.2,imag(symgray(k))+.3,...
                   %      num2str(mapgray(k)));
                  %  
                 %   text(real(symbin(k))-0.2,imag(symbin(k))-.15,...
                %        dec2base(mapbin(k),2,4),'Color',[1 0 0]);
               %     text(real(symbin(k))-0.2,imag(symbin(k))-.3,...
              %          num2str(mapbin(k)),'Color',[1 0 0]);
             %   end
            %end
        else
            error("QPSK et M différent de 4 ? problème de définition");
        end
    case 'QAM'
        symboles_envoye = qammod(info_binaire,M,'InputType','bit');
    case 'PSK'
        info_dec = bit2int(info_binaire,log2(M));
        symboles_envoye = pskmod(info_dec,M,0);
    otherwise
        error("Type non reconnu ou non supporté : doit être ASK ou QPSK ou QAM ou PSK");
end


%% Modulation
Suite_diracs = kron(symboles_envoye, [1 zeros(1, Ns-1)]);
Suite_diracs_decale=[Suite_diracs zeros(1,floor(N/2))]; 
xe_decale = filter(h, 1, Suite_diracs_decale);
xe = xe_decale(floor(N/2)+1:end);

%% Bruit
if E_bN0Db < 0
    x_bruite = xe;
else
    P_re =  mean(abs(xe).^2);
    Sigma_n = sqrt((P_re*2*Fe/Rb)/(2*log2(M)*10.^(E_bN0Db/10)));
    bruit = Sigma_n*randn(1, length(xe))+1i*Sigma_n*randn(1, length(xe));
    x_bruite = xe + bruit;
end

%% Demodulation
x_bruite_decale = [x_bruite zeros(1,floor(N/2))];
z_decale = filter(hr, 1, x_bruite_decale);
z = z_decale(floor(N/2)+1:end);

z_echant = z(n0:Ns:end);

facteur = n0;
if affichage
    scatterplot(z_echant);
end
switch type
    case 'ASK'
        symbole_recu = ((abs(z_echant/facteur)<(mapping(3)+mapping(4))/2)*mapping(3) + (abs(z_echant/facteur)>=(mapping(3)+mapping(4))/2)*mapping(4)).*sign(z_echant);
        information_sortie_reshape=[symbole_recu>(mapping(2)+mapping(3))/2 ; abs(symbole_recu)>(mapping(3)+mapping(4))/2];
        information_sortie = reshape(information_sortie_reshape,1,nb_bits);
    case {'QPSK','QAM'}
        
        information_sortie = reshape(qamdemod(z_echant,M,'OutputType','bit'), 1, nb_bits);
    case 'PSK'
       
       information_sortie = reshape(int2bit(pskdemod(z_echant,M),log2(M)), 1, nb_bits);
end
end