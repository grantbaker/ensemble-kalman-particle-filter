clear all
close all
clc

global H R gamma;

gamma = 0.7;

% Constants associated with the filter
ensNum = 30;
varNum = 3;

% Set values for true solution
deltaT = 10^-4;
tFin = 10 + deltaT; % adding extra deltaT since we don't want to observe system at time = 0
numObs = 200;
stepsBetweenObs = (tFin - deltaT)/(deltaT*numObs);
InitialCond = randn(varNum,1);

% Obtain trajoctory of system
[T, XT] = SDESolver(deltaT, 1, tFin, InitialCond);
XT = XT(2:end,:)';T = T(2:end);

XTObserved = XT(:,1:stepsBetweenObs:end);

% Initialize observation matrix and obs err cov mat
sigma0 = sqrt(1/10);
tmp = eye(varNum);
H = tmp(2:3,:);
%H = tmp;
R = (1/10)*eye(size(H,1));
% Obtain observations
obs = H*XTObserved + sigma0*randn(size(H,1),200);

% Initialize ensemble with random numbers
EnKPF = randn(varNum,ensNum);
% Store the results
EnKPF_store = zeros(varNum, ensNum, 200);

for ii = 1:200
    disp(ii)
    % EnKPF
    % forecast ensemble
    [~,sol] = SDESolver(deltaT, ensNum, stepsBetweenObs*deltaT, EnKPF');
    EnKPF_store(:,:,ii) = squeeze(sol(end,:,:))';
    
    % analysis update
    update = EnKPF_update(EnKPF_store(:,:,ii), obs(:,ii), varNum, ensNum);
                        
    EnKPF = update;
    EnKPF_store(:,:,ii) = EnKPF;
end

mean_EnKPF = squeeze(mean(EnKPF_store,2));
RMS = sqrt(mean((XTObserved-mean_EnKPF).^2,1));
RMS_unobserved = sqrt(mean((XTObserved(1,:) - mean_EnKPF(1,:)).^2,1));

figure
hold on
plot(mean_EnKPF(1,:))
plot(XTObserved(1,:))
hold off
