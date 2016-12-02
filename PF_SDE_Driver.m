clear all
close all
clc

% Constants associated with the filter
    varNum = 3;

% Initialize the particles and weights
    NP = 300; % Number of particles
    wPF = ones(1,NP)/NP;
    wPF_store = zeros(NP, 200);
    vPF = randn(varNum, NP);
    vPF_store = randn(varNum, NP, 200);

% Set values for true solution
    deltaT = 10^-4;
    % adding extra deltaT since we don't want to observe system at time = 0 
    tFin = 10 + deltaT;
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
    
for ii = 1:200
    disp(ii)
    % PF
    % forecast particles
    [~,sol] = SDESolver(deltaT, NP, stepsBetweenObs*deltaT, vPF');
    vPF = squeeze(sol(end,:,:))';
    
    for jj = 1:NP
        wPF(jj) = exp(-.5*(norm(obs(:,ii) - H*vPF(:,jj),2)/sigma0)^2);
    end
    wPF = wPF/sum(wPF);
    wPF_store(:,ii) = wPF;

    % Resample
    NN = randsample(NP,NP,true,wPF);
    vPF = vPF(:,NN);
    vPF_store(:,:,ii) = vPF;
end

% Compute the mean at each time
vPF_m = zeros(varNum, 200);
for ii = 1:varNum
   vPF_m(ii, :) = sum(squeeze(vPF_store(ii,:,:)).*wPF_store,1);
end

RMS = sqrt(mean((vPF_m - XTObserved).^2, 1));
RMS_unobserved = sqrt(mean((XTObserved(1,:) - vPF_m(1,:)).^2,1));

figure
hold on
plot(vPF_m(1,:))
plot(XTObserved(1,:))
hold off
