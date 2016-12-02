clear all
close all
clc

global H R gamma;

gamma = 0.95;

% Constants associated with the filter
    ensNum = 30;
    varNum = 3;
    deltaT = 10^-4;
    tFin = 10 + deltaT; % adding extra deltaT since we don't want to observe system at time = 0
    numObs = 200;
    numberOfSteps = (tFin - deltaT)/deltaT;
    stepsBetweenObs = (tFin - deltaT)/(deltaT*numObs);
    InitialCond = randn(varNum,1);

% Obtain trajoctory of system
    [T, XT] = SDESolver(deltaT, 1, tFin, InitialCond);
    XT = XT(2:end,:)';T = T(2:end);

    TObserved = T(:,1:stepsBetweenObs:end);
    XTObserved = XT(:,1:stepsBetweenObs:end);

% Initialize observation matrix and obs err cov mat
    sigma0 = sqrt(1/5);
    tmp = eye(varNum);
    H = tmp(2:3,:);
    %H = tmp;
    R = (1/10)*eye(size(H,1));
    % Obtain observations
    obs = H*XTObserved + sigma0*randn(size(H,1),numObs);

% Initialize ensemble
    EnKPF = bsxfun(@plus, XTObserved(:,1), sigma0*randn(varNum, ensNum));
% Store the results
    EnKPF_store = zeros(varNum, ensNum, numberOfSteps);

for ii = 1:numObs
    disp(ii)
    % EnKPF
    % forecast ensemble
    [~,sol] = SDESolver(deltaT, ensNum, stepsBetweenObs*deltaT, EnKPF');
    EnKPF_store(:,:,((ii-1)*stepsBetweenObs+1):(ii*stepsBetweenObs)) = ...
                permute(sol(:,:,:), [3,2,1]);
    EnKPF = squeeze(sol(end,:,:))';
    
    % analysis update
    update = EnKPF_update(EnKPF, obs(:,ii), varNum, ensNum);
                        
    EnKPF = update;
    EnKPF_store(:,:,ii*stepsBetweenObs) = EnKPF;
end

mean_EnKPF = squeeze(mean(EnKPF_store,2));
RMS = sqrt(mean((XT - mean_EnKPF).^2,1));
RMS_unobserved = sqrt(mean((XT(1,:) - mean_EnKPF(1,:)).^2,1));

figure
hold on
plot(T, mean_EnKPF(1,:))
plot(T, XT(1,:))
legend('EnKPF', 'True')
hold off
