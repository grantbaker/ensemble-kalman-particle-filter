clear all
close all
clc

varNum = 4;

% Initialize the particles and weights
    NP = 500; % Number of particles
    totalTimeSteps = 10000;
    numObs = 200;
    stepsBetweenObs = totalTimeSteps/numObs;
    tFin = 50;
    deltaT = tFin/totalTimeSteps;
    wPF = ones(1,NP)/NP;
    wPF_store = zeros(NP, totalTimeSteps);
    vPF = randn(varNum, NP);
    vPF_store = randn(varNum, NP, totalTimeSteps);

% Obtain trajoctory of system
    [T,XT] = ode45(@RHS_L96,linspace(0,tFin,totalTimeSteps+1),randn(varNum,1));
    XT = XT(2:end,:)';T = T(2:end);
    
    XTObserved = XT(:,1:stepsBetweenObs:end);
    
% Initialize observation matrix and obs err cov mat
    sigma0 = sqrt(1/100);
    tmp = eye(varNum);
    H = tmp(2:2:end,:);
    %H = tmp;
    R = (1/10)*eye(size(H,1));
    % Obtain observations
    obs = H*XTObserved + sigma0*randn(size(H,1),numObs);
    
for ii = 1:numObs
    disp(ii)
    % PF
    % forecast Particles
    for k = 1:NP
        [~,sol] = ode45(@RHS_L96, 0:deltaT:(stepsBetweenObs*deltaT), vPF(:,k));
        vPF_store(:,k,((ii-1)*stepsBetweenObs+1):(ii*stepsBetweenObs)) = ...
                permute(sol(1:(end-1),:), [2,1]);
        vPF(:,k) = squeeze(sol(end,:))';
    end
    
    for jj = 1:NP
        wPF(jj) = exp(-.5*(norm(obs(:,ii) - H*vPF(:,jj),2)/sigma0)^2);
    end
    wPF = wPF/sum(wPF);
    wPF_store(:,((ii-1)*stepsBetweenObs+1):(ii*stepsBetweenObs)) = ...
                    repmat(wPF,stepsBetweenObs,1)';
    % Resample
    NN = randsample(NP,NP,true,wPF);
    vPF = vPF(:,NN);
    vPF_store(:,:,ii*stepsBetweenObs) = vPF;
end

% Compute the mean at each time
vPF_m = zeros(varNum, 200);
for ii = 1:varNum
   vPF_m(ii, :) = sum(squeeze(vPF_store(ii,:,:)).*wPF_store,1);
end

RMS = sqrt(mean((vPF_m - XT).^2, 1));
