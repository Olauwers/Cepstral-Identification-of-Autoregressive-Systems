% Consider the influence of a non-white known input signal.

%% Initial test on synthetic data:
% Generate a test model:
sys = zpk([],[0.5,0.5j,-0.5j],1,1); % Third order AR model.

poles = pole(sys); % Poles: 0.5, 0.5j, and -0.5j
transf = tf(sys); % TF: a(z) = z^3 - 0.5 z^2 + 0.25 z - 0.125

%% Generate the random input/output-data:
N = 5*10^3;
freq = 1.20;
colorednoise = dsp.ColoredNoise('Color','custom', 'InverseFrequencyPower', freq,'SamplesPerFrame',N);
input = colorednoise();
output = lsim(sys,input);

%% Calculate the poles:
cepsin = ifft(log(pwelch(input,[],[],'twosided')),'symmetric');
cepsout = ifft(log(pwelch(output,[],[],'twosided')),'symmetric');
cepsmodel = cepsout-cepsin;
sysceps = tf(1,cepsarid(cepsmodel,3)',-1);
syscepswhite = tf(1,cepsarid(cepsout,3)',-1);
polesceps = pole(sysceps);
polescepswhite = pole(syscepswhite);

%% Visualize the results:
figure(1)
clf
plot(input)
title('Data')

figure(2)
clf
hold on
pzmap(sys)
plot(real(polesceps),imag(polesceps),'*')
plot(real(polescepswhite),imag(polescepswhite),'*')
hold off
title('Poles')
legend('original','cepstrum','cepstrum white')

[mag1,~,wout1] = bode(sys);
[mag2,~,wout2] = bode(sysceps);
[mag3,~,wout3] = bode(syscepswhite);

figure(3)
clf
hold on
plot(wout1,20*log10(squeeze(mag1)))
plot(wout2,20*log10(squeeze(mag2)))
plot(wout3,20*log10(squeeze(mag3)))
hold off
title('Bode diagram')
legend('original','cepstrum','cepstrum white')

F = linspace(0,0.5,1000);
PSD = 1./abs((F.^freq));
figure(4)
clf
plot(2*pi*F(2:end),mag2db(PSD(2:end)))
title('Power spectral density')
set(gca, 'XScale', 'log')