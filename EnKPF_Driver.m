clear all
close all
clc

global H R gamma;

gamma = 0.7;

% Constants associated with the filter
ensNum = 30;
varNum = 5;

% Obtain trajoctory of system
[T,XT] = ode45(@RHS_L96,linspace(0,50,201),randn(varNum,1));
XT = XT(2:end,:)';T = T(2:end);

% Initialize observation matrix and obs err cov mat
sigma0 = sqrt(1/10);
tmp = eye(varNum);
H = tmp(2:2:end,:);
%H = tmp;
R = (1/10)*eye(size(H,1));
% Obtain observations
obs = H*XT + sigma0*randn(size(H,1),200);

% Initialize ensemble with random numbers
EnKPF = randn(varNum,ensNum);
% Store the results
EnKPF_store = zeros(varNum, ensNum, 200);

for ii = 1:200
    % EnKPF
    % forecast ensemble
    for k = 1:ensNum
        [~,sol] = ode45(@RHS_L96, [0 0.25], EnKPF(:,k));
        EnKPF_store(:,k,ii) = squeeze(sol(end,:))';
    end
    
    % analysis update
    update = EnKPF_update(EnKPF_store(:,:,ii), obs(:,ii), varNum, ensNum);
                        
    EnKPF = update;
    EnKPF_store(:,:,ii) = EnKPF;
end

mean_EnKPF = squeeze(mean(EnKPF_store,2));
RMS = sqrt(mean((XT-mean_EnKPF).^2,1));