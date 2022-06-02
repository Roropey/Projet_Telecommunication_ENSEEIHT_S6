function [info_entree,info_recu,x,y_bruite,z,z_echant] = Propagation_Multi_Canal(info_entree,nb_bits,Fe,Rb,n0,a,h,hr,alpha_0,tau_0,alpha_1,tau_1,E_bN0dB)




%% Variables initiales
if (info_entree == Inf)
    info_entree = randi([0,1], 1,nb_bits);
end;

%% Modulateur

% Variables
Ns = log2(length(a))*Fe/Rb;
Ts = Ns/Fe;

% Calculs
mapping = info_entree.*(a(2) - a(1)) + a(1);
Suite_diracs = kron(mapping, [1 zeros(1, Ns-1)]);
x = filter(h, 1, Suite_diracs);


%% Canal

hc = zeros(1,(max(tau_0,tau_1)+Ts)*Fe);
hc(tau_0*Fe+1) = alpha_0;
hc(tau_1*Fe+1) = alpha_1;

y=filter(hc,1,x);

%% Bruit
if (E_bN0dB==Inf)
    y_bruite = y;
else
    P_re =  mean(abs(y).^2);
    Sigma_n = sqrt((P_re*Ns)/(2*log2(length(a))*10.^(E_bN0dB/10)));
    bruit = Sigma_n*randn(1, length(y))+1i*Sigma_n*randn(1, length(y));
    y_bruite = y + bruit;
end

%% Demodulateur

%Filtrage de réception
z = filter(hr, 1, y_bruite);

z_echant = z(n0:Ns:end);
info_recu = real(z_echant) > 0;

end