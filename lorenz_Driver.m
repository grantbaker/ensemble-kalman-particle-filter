clear all
close all
clc

global H R;

% Constants associated with the filter
ensNum = 3;
varNum = 3;
numCenters = 2;
N = 1;

% Obtain trajoctory of system
[T,XT] = ode45(@lorenz,linspace(0,50,201),randn(varNum,1));
XT = XT(2:end,:)';T = T(2:end);

% Initialize observation matrix and obs err cov mat
sigma0 = sqrt(1/10);
tmp = eye(varNum);
H = tmp(1:2:end,:);
R = (1/10)*eye(size(H,1));
% Obtain observations
obs = H*XT + sigma0^2*randn(size(H,1),200);

% Initialize mixing probability randomly
mixingProb = unifrnd(0, 1, [numCenters,1]);
mixingProb = mixingProb/sum(mixingProb);

% Initialize ensemble with random numbers
XEnKF = randn(varNum,ensNum);
% Store the results
XEnKF_store = zeros(varNum, ensNum, 200);

for ii = 1:200
    % XEnKF
    % forecast ensemble
    for k = 1:ensNum
        [~,sol] = ode45(@lorenz, [0 0.25], XEnKF(:,k));
        XEnKF_store(:,k,ii) = squeeze(sol(end,:))';
    end
    
    % analysis update
    [update, mixingProb] = XEnKF_update(XEnKF_store(:,:,ii), obs(:,ii), ...
                            mixingProb, varNum, numCenters, N, ensNum);
                        
    XEnKF = update;
    XEnKF_store(:,:,ii) = XEnKF;
end

mean_XEnKF = squeeze(mean(XEnKF_store,2));
RMS = sqrt(mean((XT-mean_XEnKF).^2,1));

figure
hold on
plot(RMS)
hold off