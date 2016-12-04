clear all
close all
clc

% Constants associated with the filter
    ensNum = 100;
    varNum = 5;
    totalTimeSteps = 5000;
    numObs = 500;
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
    H = tmp(2:2:end,:);
    %H = tmp;
    R = (1/10)*eye(size(H,1));
    % Obtain observations
    obs = H*XTObserved + sigma0*randn(size(H,1),numObs);

% Initialize ensemble with random numbers
    EnKF = randn(varNum,ensNum);
    % Store the results
    EnKF_store = zeros(varNum, ensNum, totalTimeSteps);

for ii = 1:numObs
    disp(ii)
    % EnKPF
    % forecast ensemble
    for k = 1:ensNum
        [~,sol] = ode45(@RHS_L96, 0:deltaT:(stepsBetweenObs*deltaT), EnKF(:,k));
        EnKF_store(:,k,((ii-1)*stepsBetweenObs+1):(ii*stepsBetweenObs)) = ...
                permute(sol(1:(end-1),:), [2,1]);
        EnKF(:,k) = squeeze(sol(end,:))';
    end
    
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
RMS_observed = sqrt(mean((XT(2:2:end,:) - mean_EnKF(2:2:end,:)).^2,1));
RMS_unobserved = sqrt(mean((XT(1:2:end,:) - mean_EnKF(1:2:end,:)).^2,1));

figure
hold on
plot(RMS)
plot(RMS_observed)
plot(RMS_unobserved)
legend('RMS of All Variables', 'RMS of Observed', 'RMS of Unobserved')
hold off