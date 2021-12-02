% Apply the cepstral system identification approach to the historical
% Yule's sunspot problem.

%% Initialize the workspace:
close all
clearvars
load('sunspot.mat')

%% Determine the coefficients of the AR model:
p = 2; % AR model order.

% Generate cepstrum coefficients:
% w = detrend(w); % Because the papers also detrend the data!
ceps = ifft(log(pmtm(w,[],[],'twosided')),'symmetric');

% Identify the models via cepstrum and least-squares:
syscepsno = tf(1,cepsarid(ceps,p)',-1);
syslsno = tf(1,[1; getpvec(ar(w,p,'ls'))]',-1);
sysburgno = tf(1,[1; getpvec(ar(w,p,'burg'))]',-1);
syslpcno = tf(1,lpc(w,p),-1);

% Add the model from the paper:
sysyule = tf(1,[1 -1.34254 0.65504],-1); % G. U. Yule 1927
v = w;
w = detrend(w); % Because the papers also detrend the data!
ceps = ifft(log(pmtm(w,[],[],'twosided')),'symmetric');

% Identify the models via cepstrum and least-squares (after detrending):
sysceps = tf(1,cepsarid(ceps,p)',-1);
sysls = tf(1,[1; getpvec(ar(w,p,'ls'))]',-1);
sysburg = tf(1,[1; getpvec(ar(w,p,'burg'))]',-1);
syslpc = tf(1,lpc(w,p),-1);
cepsarid(ceps,p)
getpvec(ar(w,p,'ls'))
getpvec(ar(w,p,'burg'))
lpc(w,p)

%% Visualize the results:
[magceps,~,woutceps] = bode(sysceps);
[magls,~,woutls] = bode(sysls);
[magburg,~,woutburg] = bode(sysburg);
[magyule,~,woutyule] = bode(sysyule);
[maglpc,~,woutlpc] = bode(syslpc);

figure(1)
clf
hold on
plot(woutyule,20*log10(squeeze(magyule)))
plot(woutceps,20*log10(squeeze(magceps)))
plot(woutls,20*log10(squeeze(magls)))
plot(woutburg,20*log10(squeeze(magburg)))
plot(woutlpc,20*log10(squeeze(maglpc)))
hold off
title('Bode diagram')
legend('Yule', 'cepstrum','leastsq','burg', 'lpc')

[magcepsno,~,woutcepsno] = bode(syscepsno);
[maglsno,~,woutlsno] = bode(syslsno);
[magburgno,~,woutburgno] = bode(sysburgno);
[maglpcno,~,woutlpcno] = bode(syslpcno);

figure(2)
clf
hold on
plot(woutyule,20*log10(squeeze(magyule)))
plot(woutcepsno,20*log10(squeeze(magcepsno)))
plot(woutlsno,20*log10(squeeze(maglsno)))
plot(woutburgno,20*log10(squeeze(magburgno)))
plot(woutlpcno,20*log10(squeeze(maglpcno)))
hold off
title('Bode diagram')
legend('Yule', 'cepstrum','leastsq','burg', 'lpc')

figure(3)
clf
hold on
plot(v)
plot(w)
plot(v-w)
hold off
title('Sunspot numbers') 
legend('original','detrended','trend')
