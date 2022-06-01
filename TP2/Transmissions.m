%% Nettoyage
close all;
clear;

%% Variables Initiales
alpha  = 0.5;
Fp = 2000;
Fe = 480000;
Rb = 48000;
nb_bits = 120000;
N = 201;
seuil_erreur = 1000;
n0 = 1;
E_bN0Db = -1;

E_bN0db = 0:0.1:6;

%(Fe,Rb,N,type,M,nb_bits,n0,h,hr,E_bN0db,seuil_erreur,TEB_th)

%% 4-ASK

TEB_th_4ASK = (2/log2(4)) * (1-1/4) * qfunc(sqrt(((6*log2(4))/(4*4-1)).*10.^(E_bN0db/10)));
[TEB_4ASK signal_emis_4ASK] = EtudeTransmission(Fe,Rb,N,'ASK',4,nb_bits,n0,alpha,E_bN0db,seuil_erreur,TEB_th_4ASK);

pause;
close all;
%% QPSK

TEB_th_QPSK = (4/ log2(4)).*(1-(1/sqrt(4))).*qfunc(sqrt(((3*log2(4))/(4-1)).*10.^(E_bN0db/10)));
[TEB_QPSK signal_emis_QPSK] = EtudeTransmission(Fe,Rb,N,'QPSK',4,nb_bits,n0,alpha,E_bN0db,seuil_erreur,TEB_th_QPSK);

pause;
close all;
%% 8-PSK

TEB_th_8PSK = (2/log2(8)) * qfunc(sqrt(2*(10.^(E_bN0db/10))*log2(8))*sin(pi/8));
[TEB_8PSK signal_emis_8PSK] = EtudeTransmission(Fe,Rb,N,'PSK',8,nb_bits,n0,alpha,E_bN0db,seuil_erreur,TEB_th_8PSK);

pause;
close all;
%% 16-QAM

TEB_th_16QAM = (4/ log2(16)).*(1-(1/sqrt(16))).*qfunc(sqrt(((3*log2(16))/(16-1)).*10.^(E_bN0db/10)));
[TEB_16QAM signal_emis_16QAM] = EtudeTransmission(Fe,Rb,N,'QAM',16,nb_bits,n0,alpha,E_bN0db,seuil_erreur,TEB_th_16QAM);

pause;
close all;
%% Comparaison

%TEB

figure('Name','Comparaison TEB');
s1_TEB = semilogy(E_bN0db,TEB_4ASK);
hold on;
s2_TEB = semilogy(E_bN0db,TEB_QPSK);
s3_TEB = semilogy(E_bN0db,TEB_8PSK);
s4_TEB = semilogy(E_bN0db,TEB_16QAM);
hold off;

xlabel('Eb/N0 (dB)');
ylabel('TEB');
title('TEBs pour différentes modulations');
legend([s1_TEB s2_TEB s3_TEB s4_TEB],"4-ASK","QPSK","8-PSK","16-QAM");

figure('Name','Comparaison TEB théorique');
s1 = semilogy(E_bN0db,TEB_th_4ASK);
hold on;
s2 = semilogy(E_bN0db,TEB_th_QPSK);
s3 = semilogy(E_bN0db,TEB_th_8PSK);
s4 = semilogy(E_bN0db,TEB_th_16QAM);
hold off;

xlabel('Eb/N0 (dB)');
ylabel('TEB');
title('TEBs pour différentes modulations');
legend([s1 s2 s3 s4],"4-ASK","QPSK","8-PSK","16-QAM");

%DSP

DSP_4ASK = fftshift(abs(fft(xcorr(signal_emis_4ASK,'unbiased'),10000)));
plage_4ASK=(-Fe/2 : Fe/(length(DSP_4ASK)-1) : Fe/2);

DSP_QPSK = fftshift(abs(fft(xcorr(signal_emis_QPSK,'unbiased'),10000)));
plage_QPSK=(-Fe/2 : Fe/(length(DSP_QPSK)-1) : Fe/2);

DSP_8PSK = fftshift(abs(fft(xcorr(signal_emis_8PSK,'unbiased'),10000)));
plage_8PSK=(-Fe/2 : Fe/(length(DSP_8PSK)-1) : Fe/2);

DSP_16QAM = fftshift(abs(fft(xcorr(signal_emis_16QAM,'unbiased'),10000)));
plage_16QAM=(-Fe/2 : Fe/(length(DSP_16QAM)-1) : Fe/2);

figure('Name','Comparaison DSP');
s1_DSP = semilogy(plage_4ASK,DSP_4ASK);
hold on;
s2_DSP = semilogy(plage_QPSK,DSP_QPSK);
s3_DSP = semilogy(plage_8PSK,DSP_8PSK);
s4_DSP = semilogy(plage_16QAM,DSP_16QAM);
hold off;

xlabel('Hz');
ylabel('DSP');
title('Densité spectrales de puissances pour différentes modulations');
legend([s1_DSP s2_DSP s3_DSP s4_DSP],"4-ASK","QPSK","8-PSK","16-QAM");
