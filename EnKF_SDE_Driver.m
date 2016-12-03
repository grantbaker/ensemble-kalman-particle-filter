clear all
close all
clc

% Constants associated with the filter
ensNum = 20;
varNum = 3;

% Set values for true solution
    deltaT = 10^-4;
    tFin = 10 + deltaT; % adding extra deltaT since we don't want to observe system at time = 0
    numObs = 2000;
    numberOfSteps = (tFin - deltaT)/deltaT;
    stepsBetweenObs = (tFin - deltaT)/(deltaT*numObs);
    InitialCond = randn(varNum,1);

% Obtain trajoctory of system
    [T, XT] = SDESolver(deltaT, 1, tFin, InitialCond);
    XT = XT(2:end,:)';T = T(2:end);

    XTObserved = XT(:,1:stepsBetweenObs:end);

% Initialize observation matrix and obs err cov mat
    sigma0 = sqrt(1/5);
    tmp = eye(varNum);
    H = tmp(2:3,:);
    %H = tmp;
    R = sigma0^2*eye(size(H,1));
    % Obtain observations
    obs = H*XTObserved + sigma0*randn(size(H,1), numObs);

% Store the results
    EnKF_store = zeros(varNum, ensNum, numberOfSteps);
    EnKF = bsxfun(@plus, XTObserved(:,1), randn(varNum, ensNum));

for ii = 1:numObs
    disp(ii)
    % EnKPF
    % forecast ensemble
    [~,sol] = SDESolver(deltaT, ensNum, stepsBetweenObs*deltaT, EnKF');
    EnKF_store(:,:,((ii-1)*stepsBetweenObs+1):(ii*stepsBetweenObs)) = ...
                permute(sol(:,:,:), [3,2,1]);
    EnKF = squeeze(sol(end,:,:))';
    
    % compute the ensemble mean
    mu = (1/ensNum)*sum(EnKF, 2);
    
    % compute the ensemble covariance
    A = (bsxfun(@plus, EnKF, - mu))/(ensNum-1);
    
    % compute the Kalman gain matrix
    K = (A*(H*A)')/((H*A*(H*A)') + R);
    
    % apply the Kalman filter update on each ensemble member
    for jj=1:ensNum
        % Perturbation for y
        eps_y = normrnd(0,sigma0);
        EnKF(:,jj) = EnKF(:,jj) + K*(obs(:,ii) + eps_y - ...
                                H*EnKF(:,jj));
    end
    
    EnKF_store(:,:,ii*stepsBetweenObs) = EnKF;
end

mean_EnKF = squeeze(mean(EnKF_store,2));
RMS = sqrt(mean((XT - mean_EnKF).^2,1));
RMS_unobserved = sqrt(mean((XT(1,:) - mean_EnKF(1,:)).^2,1));

figure
hold on
plot(T, mean_EnKF(1,:))
plot(T, XT(1,:))
legend('EnKF', 'True')
hold off
