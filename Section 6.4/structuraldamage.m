% Apply the cepstral system identification approach to the a practical
% identification problem from structural health monitoring.

%% Initialize the workspace:
% close all
clearvars
load('structuraldamage.mat')
set(0,'DefaultFigureWindowStyle','normal')

p = 3; % Order of the AR-model.
N = 50; % Size of the window.

number = floor(length(data{1})/N);
used = number*N;

%% Via the cepstral system identification approach:
% Measurements of undamaged structure (1):
coef1 = zeros(number,p);
measurement1 = data{1};
for i = 1:number
    ceps = ifft(log(pmtm(measurement1((i-1)*N+1:i*N),[],[],'twosided')),'symmetric');
    vec = cepsarid(ceps,p)';
    coef1(i,:) = vec(2:end);
end

% Measurements of undamaged structure (2):
coef2 = zeros(number,p);
measurement2 = data{2};
for i = 1:number
    ceps = ifft(log(pmtm(measurement2((i-1)*N+1:i*N),[],[],'twosided')),'symmetric');
    vec = cepsarid(ceps,p)';
    coef2(i,:) = vec(2:end);
end

% Measurements of damaged structure (3):
coef3 = zeros(number,p);
measurement3 = data{3};
for i = 1:number
    ceps = ifft(log(pmtm(measurement3((i-1)*N+1:i*N),[],[],'twosided')),'symmetric');
    vec = cepsarid(ceps,p)';
    coef3(i,:) = vec(2:end);
end

% Visualization:
figure(1)
clf
hold on
plot(coef1(:,1))
plot(coef2(:,1))
plot(coef3(:,1))
hold off

figure(2)
clf
hold on
plot(coef1(:,2))
plot(coef2(:,2))
plot(coef3(:,2))
hold off

figure(3)
clf
hold on
plot(coef1(:,3))
plot(coef2(:,3))
plot(coef3(:,3))
hold off

%% Via the least-squares system identification approach:
% Measurements of undamaged structure (1):
coef1 = zeros(number,p);
measurement1 = data{1};
for i = 1:number
    vec = getpvec(ar(measurement1((i-1)*N+1:i*N),p,'ls'))';
    coef1(i,:) = vec(1:end);
end

% Measurements of undamaged structure (2):
coef2 = zeros(number,p);
measurement2 = data{2};
for i = 1:number
    vec = getpvec(ar(measurement2((i-1)*N+1:i*N),p,'ls'))';
    coef2(i,:) = vec(1:end);
end

% Measurements of damaged structure (3):
coef3 = zeros(number,p);
measurement3 = data{3};
for i = 1:number
    vec = getpvec(ar(measurement3((i-1)*N+1:i*N),p,'ls'))';
    coef3(i,:) = vec(1:end);
end

% Visualization:
figure(4)
clf
hold on
plot(coef1(:,1))
plot(coef2(:,1))
plot(coef3(:,1))
hold off

figure(5)
clf
hold on
plot(coef1(:,2))
plot(coef2(:,2))
plot(coef3(:,2))
hold off

figure(6)
clf
hold on
plot(coef1(:,3))
plot(coef2(:,3))
plot(coef3(:,3))
hold off

%% Via the Burg system identification approach:
% Measurements of undamaged structure (1):
coef1 = zeros(number,p);
measurement1 = data{1};
for i = 1:number
    vec = getpvec(ar(measurement1((i-1)*N+1:i*N),p,'burg'))';
    coef1(i,:) = vec(1:end);
end

% Measurements of undamaged structure (2):
coef2 = zeros(number,p);
measurement2 = data{2};
for i = 1:number
    vec = getpvec(ar(measurement2((i-1)*N+1:i*N),p,'burg'))';
    coef2(i,:) = vec(1:end);
end

% Measurements of damaged structure (3):
coef3 = zeros(number,p);
measurement3 = data{3};
for i = 1:number
    vec = getpvec(ar(measurement3((i-1)*N+1:i*N),p,'burg'))';
    coef3(i,:) = vec(1:end);
end

% Visualization:
figure(7)
clf
hold on
plot(coef1(:,1))
plot(coef2(:,1))
plot(coef3(:,1))
hold off

figure(8)
clf
hold on
plot(coef1(:,2))
plot(coef2(:,2))
plot(coef3(:,2))
hold off

figure(9)
clf
hold on
plot(coef1(:,3))
plot(coef2(:,3))
plot(coef3(:,3))
hold off
 
%% Via the LPC system identification approach:
% Measurements of undamaged structure (1):
coef1 = zeros(number,p);
measurement1 = data{1};
for i = 1:number
    vec = lpc(measurement1((i-1)*N+1:i*N),p);
    coef1(i,:) = vec(2:end);
end

% Measurements of undamaged structure (2):
coef2 = zeros(number,p);
measurement2 = data{2};
for i = 1:number
    vec = lpc(measurement2((i-1)*N+1:i*N),p);
    coef2(i,:) = vec(2:end);
end

% Measurements of damaged structure (3):
coef3 = zeros(number,p);
measurement3 = data{3};
for i = 1:number
    vec = lpc(measurement3((i-1)*N+1:i*N),p);
    coef3(i,:) = vec(2:end);
end

% Visualization:
figure(10)
clf
hold on
plot(coef1(:,1))
plot(coef2(:,1))
plot(coef3(:,1))
hold off

figure(11)
clf
hold on
plot(coef1(:,2))
plot(coef2(:,2))
plot(coef3(:,2))
hold off

figure(12)
clf
hold on
plot(coef1(:,3))
plot(coef2(:,3))
plot(coef3(:,3))
hold off