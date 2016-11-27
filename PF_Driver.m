clear all
close all
clc

varNum = 5;

% Obtain trajoctory of system
    [T,XT] = ode45(@RHS_L96,linspace(0,50,201),randn(varNum,1));
    XT = XT(2:end,:)';T = T(2:end);

% Initialize the particles and weights
    NP = 100; % Number of particles
    wPF = ones(1,NP)/NP;
    wPF_store = zeros(NP, 200);
    vPF = randn(varNum, NP);
    vPF_store = randn(varNum, NP, 200);
    
% Initialize observation matrix and obs err cov mat
    sigma0 = sqrt(1/100);
    tmp = eye(varNum);
    H = tmp(2:2:end,:);
    %H = tmp;
    R = (1/10)*eye(size(H,1));
    % Obtain observations
    obs = H*XT + sigma0*randn(size(H,1),200);
    
for ii = 1:200
    % PF
    % forecast Particles
    for k = 1:NP
        [~,sol] = ode45(@RHS_L96, [0 0.25], vPF(:,k));
        vPF(:,k) = squeeze(sol(end,:))';
    end
    
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
RMS = sqrt(mean((vPF_m - XT).^2, 1));