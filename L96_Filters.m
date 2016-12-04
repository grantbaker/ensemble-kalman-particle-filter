clear all
close all
clc

global H R gamma;

gamma = 0.7;

% Constants associated with the filter
    ensNum = 100;
    varNum = 5;
    totalTimeSteps = 5000;
    numObs = 200; %Needs to evenly divide totalTimeSteps
    stepsBetweenObs = totalTimeSteps/numObs;
    tFin = 25;
    deltaT = tFin/totalTimeSteps;

% Obtain trajoctory of system
    [T,XT] = ode45(@RHS_L96,linspace(0,tFin,totalTimeSteps+1),randn(varNum,1));
    XT = XT(2:end,:)';T = T(2:end);
    
    XTObserved = XT(:,1:stepsBetweenObs:end);

% Initialize observation matrix and obs err cov mat
    sigma0 = sqrt(1/10);
    tmp = eye(varNum);
    H = tmp(1:3:end,:);
    %H = tmp;
    R = (1/10)*eye(size(H,1));
    % Obtain observations
    obs = H*XTObserved + sigma0*randn(size(H,1),numObs);

% Initialize ensemble for EnKPF with random numbers
    EnKPF = randn(varNum,ensNum);
    % Store the results
    EnKPF_store = zeros(varNum, ensNum, totalTimeSteps);
    
% Initialize ensemble for EnKF
    EnKF = randn(varNum,ensNum);
    % Store the results
    EnKF_store = zeros(varNum, ensNum, totalTimeSteps);
    

for ii = 1:numObs
    disp(ii)
    % EnKPF
    % forecast ensemble
    for k = 1:ensNum
        [~,sol_EnKPF] = ode45(@RHS_L96, 0:deltaT:(stepsBetweenObs*deltaT), EnKPF(:,k));
        [~,sol_EnKF] = ode45(@RHS_L96, 0:deltaT:(stepsBetweenObs*deltaT), EnKF(:,k));
        
        EnKPF_store(:,k,((ii-1)*stepsBetweenObs+1):(ii*stepsBetweenObs)) = ...
                permute(sol_EnKPF(1:(end-1),:), [2,1]);
        EnKPF(:,k) = squeeze(sol_EnKPF(end,:))';
        
        EnKF_store(:,k,((ii-1)*stepsBetweenObs+1):(ii*stepsBetweenObs)) = ...
                permute(sol_EnKF(1:(end-1),:), [2,1]);
        EnKF(:,k) = squeeze(sol_EnKF(end,:))';
    end
    
    % EnKPF update
    update = EnKPF_update(EnKPF, obs(:,ii), varNum, ensNum);
                        
    EnKPF = update;
    EnKPF_store(:,:,ii*stepsBetweenObs) = EnKPF;
    
    % EnKF Update
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
RMS_EnKF = sqrt(mean((XT - mean_EnKF).^2,1));

mean_EnKPF = squeeze(mean(EnKPF_store,2));
RMS_EnKPF = sqrt(mean((XT - mean_EnKPF).^2,1));

figure
hold on
plot(T, RMS_EnKF)
plot(T, RMS_EnKPF)
legend('RMS EnKF', 'RMS EnKPF')
hold off
