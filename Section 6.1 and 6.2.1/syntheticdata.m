% Use the cepstral system identification approach on synthetic data.

%% Initialization of the workspace:
close all
clearvars
runs = 3;

%% Initial test on synthetic data:
% Generate a test model:
sys = zpk([],[0.5,0.5j,-0.5j],1,1); % Third order AR model.

poles = pole(sys); % Poles: 0.5, 0.5j, and -0.5j
transf = tf(sys); % TF: a(z) = z^3 - 0.5 z^2 + 0.25 z - 0.125

% Generate the random input/output-data:
N = 5*10^5;
t = 0:1:N-1;
input = randn(size(t));
output = lsim(sys,input);

% Calculate the poles:
ceps = ifft(log(pwelch(output,[],[],'twosided')),'symmetric');
sysceps = tf(1,cepsarid(ceps,3)',-1);
polesceps = pole(sysceps);

% Visualize the results:
figure(1)
clf
hold on
pzmap(sys)
plot(real(polesceps),imag(polesceps),'*')
hold off
title('Poles')
legend('original','cepstrum')

[mag1,~,wout1] = bode(sys);
[mag2,~,wout2] = bode(sysceps);

figure(2)
clf
hold on
plot(wout1,20*log10(squeeze(mag1)))
plot(wout2,20*log10(squeeze(mag2)))
hold off
title('Bode diagram')
legend('original','cepstrum')

%% Convergence test on synthetic data:
% Compare for different data lengths:
datarange = 2.^(4:1:20);
eceps = zeros(length(datarange),runs);
els = zeros(length(datarange),runs);
eburg = zeros(length(datarange),runs);
elpc = zeros(length(datarange),runs);
for k = 1:length(datarange)
    N = datarange(k);
    t = 0:1:N-1;
    for l = 1:runs
        input = randn(size(t));
        output = lsim(sys,input);
        ceps = ifft(log(pwelch(output,[],[],'twosided')),'symmetric');
        eceps(k,l) = norm(tf(1,cepsarid(ceps,3)',-1) - sys);
        els(k,l) = norm(tf(1,[1; getpvec(ar(output,3,'ls'))]',-1) - sys);
        eburg(k,l) = norm(tf(1,[1; getpvec(ar(output,3,'burg'))]',-1) - sys);
        elpc(k,l) = norm(tf(1,lpc(output,3),-1) - sys);
    end
end
stdceps = std(eceps,0,2);
meanceps = mean(eceps,2);
stdls = std(els,0,2);
meanls = mean(els,2);
stdburg = std(eburg,0,2);
meanburg = mean(eburg,2);
stdlpc = std(elpc,0,2);
meanlpc = mean(elpc,2);

% Visualize the results:
figure(3)
clf
semilogx(datarange,[meanceps meanls meanburg meanlpc meanceps+stdceps meanceps-stdceps meanls+stdls meanls-stdls meanburg+stdburg meanburg-stdburg meanlpc+stdlpc meanlpc-stdlpc])
title('Convergence of all algorithms')
legend('cepstrum','least','burg','lpc')
xlabel('N')
ylabel('H2-norm of the error system')

%% Noise test on synthetic data:
% Generate the random input/output-data:
N = 2^20;
t = 0:1:N-1;
input = randn(size(t));
output = lsim(sys,input);

% Calculate the poles for different SNRs:
output1 = output + 0*randn(size(t))'; % Inf
ceps1 = ifft(log(pwelch(output1,[],[],'twosided')),'symmetric');
sysceps1 = tf(1,cepsarid(ceps1,3)',-1);

output2 = output + 0.01*randn(size(t))'; % 40 dB
ceps2 = ifft(log(pwelch(output2,[],[],'twosided')),'symmetric');
sysceps2 = tf(1,cepsarid(ceps2,3)',-1);

output3 = output + 0.1*randn(size(t))'; % 20 dB
ceps3 = ifft(log(pwelch(output3,[],[],'twosided')),'symmetric');
sysceps3 = tf(1,cepsarid(ceps3,3)',-1);

output4 = output + 0.30*randn(size(t))'; % 10.4576 dB
ceps4 = ifft(log(pwelch(output4,[],[],'twosided')),'symmetric');
sysceps4 = tf(1,cepsarid(ceps4,3)',-1);

output5 = output + 0.5*randn(size(t))'; % 6.0206 dB
ceps5 = ifft(log(pwelch(output5,[],[],'twosided')),'symmetric');
sysceps5 = tf(1,cepsarid(ceps5,3)',-1);

output6 = output + 0.891*randn(size(t))'; % 1.0024 dB
ceps6 = ifft(log(pwelch(output6,[],[],'twosided')),'symmetric');
sysceps6 = tf(1,cepsarid(ceps6,3)',-1);

% Visualize the results:
[mag1,~,wout1] = bode(sys);
[mag2,~,wout2] = bode(sysceps1);
[mag3,~,wout3] = bode(sysceps2);
[mag4,~,wout4] = bode(sysceps3);
[mag5,~,wout5] = bode(sysceps4);
[mag6,~,wout6] = bode(sysceps5);
[mag7,~,wout7] = bode(sysceps6);
figure(4)
clf
hold on
plot(wout1,20*log10(squeeze(mag1)))
plot(wout2,20*log10(squeeze(mag2)))
plot(wout3,20*log10(squeeze(mag3)))
plot(wout4,20*log10(squeeze(mag4)))
plot(wout5,20*log10(squeeze(mag5)))
plot(wout6,20*log10(squeeze(mag6)))
plot(wout7,20*log10(squeeze(mag7)))
hold off
title('Bode diagram')
legend('original','Inf','40dB','20dB','10db','6dB','1dB')

%% Test influence of wrong model order:
% Generate a test model:
sysupdated = zpk([],[0.5,0.5j,-0.5j,-0.25+0.25j,-0.25-0.25j],1,1); % Fifth order AR model.

polesupdated = pole(sysupdated); % Poles: 0.5, 0.5j, -0.5j, -0.25 + 0.25j, and -0.25 - 0.25j.
transfupdated = tf(sysupdated); % TF: a(z) = z^3 - 0.5 z^2 + 0.25 z - 0.125

% Generate the random input/output-data:
N = 5*10^5;
t = 0:1:N-1;
input = randn(size(t));
output = lsim(sysupdated,input);
ceps = ifft(log(pwelch(output,[],[],'twosided')),'symmetric');

% Calculate the poles (correct model order):
syscorrect = tf(1,cepsarid(ceps,5)',-1);
polescorrect = pole(syscorrect);

% Calculate the poles (wrong model order):
syswrong1 = tf(1,cepsarid(ceps,4)',-1);
poleswrong1 = pole(syswrong1);

% Calculate the poles (wrong model order):
syswrong2 = tf(1,cepsarid(ceps,3)',-1);
poleswrong2 = pole(syswrong2);

% Calculate the poles (wrong model order):
syswrong3 = tf(1,cepsarid(ceps,6)',-1);
poleswrong3 = pole(syswrong3);

% Calculate the poles (wrong model order):
syswrong4 = tf(1,cepsarid(ceps,7)',-1);
poleswrong4 = pole(syswrong4);

% Visualize the results:
figure(5)
clf
hold on
pzmap(sysupdated)
plot(real(poleswrong4),imag(poleswrong4),'*')
plot(real(poleswrong3),imag(poleswrong3),'*')
plot(real(polescorrect),imag(polescorrect),'*')
plot(real(poleswrong1),imag(poleswrong1),'*')
plot(real(poleswrong2),imag(poleswrong2),'*')
hold off
title('Poles')
legend('original', 'seventh order', 'sixth order', 'fifth order', 'fourth order', 'third order')

[mag1,~,wout1] = bode(sysupdated);
[mag2,~,wout2] = bode(syscorrect);
[mag3,~,wout3] = bode(syswrong1);
[mag4,~,wout4] = bode(syswrong2);
[mag5,~,wout5] = bode(syswrong3);
[mag6,~,wout6] = bode(syswrong4);

figure(6)
clf
hold on
plot(wout6,20*log10(squeeze(mag6)))
plot(wout5,20*log10(squeeze(mag5)))
plot(wout1,20*log10(squeeze(mag1)))
plot(wout2,20*log10(squeeze(mag2)))
plot(wout3,20*log10(squeeze(mag3)))
plot(wout4,20*log10(squeeze(mag4)))
hold off
title('Bode diagram')
legend('original', 'seventh order', 'sixth order', 'fifth order', 'fourth order', 'third order')

% Compare for different data lengths:
datarange = 2.^((8:2:40)./2);
eceps7 = zeros(length(datarange),runs);
eceps6 = zeros(length(datarange),runs);
eceps5 = zeros(length(datarange),runs);
eceps4 = zeros(length(datarange),runs);
eceps3 = zeros(length(datarange),runs);
for k = 1:length(datarange)
    N = datarange(k);
    t = 0:1:N-1;
    for l = 1:runs
        input = randn(size(t));
        output = lsim(sysupdated,input);
        ceps = ifft(log(pwelch(output,[],[],'twosided')),'symmetric');
        eceps7(k,l) = norm(tf(1,cepsarid(ceps,7)',-1) - sysupdated);
        eceps6(k,l) = norm(tf(1,cepsarid(ceps,6)',-1) - sysupdated);
        eceps5(k,l) = norm(tf(1,cepsarid(ceps,5)',-1) - sysupdated);
        eceps4(k,l) = norm(tf(1,cepsarid(ceps,4)',-1) - sysupdated);
        eceps3(k,l) = norm(tf(1,cepsarid(ceps,3)',-1) - sysupdated);
    end
end
stdceps7 = std(eceps7,0,2);
meanceps7 = mean(eceps7,2);
stdceps6 = std(eceps6,0,2);
meanceps6 = mean(eceps6,2);
stdceps5 = std(eceps5,0,2);
meanceps5 = mean(eceps5,2);
stdceps4 = std(eceps4,0,2);
meanceps4 = mean(eceps4,2);
stdceps3 = std(eceps3,0,2);
meanceps3 = mean(eceps3,2);

% Visualize the results:
figure(7)
clf
semilogx(datarange,[meanceps7 meanceps6 meanceps5 meanceps4 meanceps3 meanceps7+stdceps7 meanceps7-stdceps7 meanceps6+stdceps6 meanceps6-stdceps6 meanceps5+stdceps5 meanceps5-stdceps5 meanceps4+stdceps4 meanceps4-stdceps4 meanceps3+stdceps3 meanceps3-stdceps3])
title('Convergence of different model orders')
legend('seventh order', 'sixth order', 'fifth order','fourth order','third order')
xlabel('N')
ylabel('H2-norm of the error system')

