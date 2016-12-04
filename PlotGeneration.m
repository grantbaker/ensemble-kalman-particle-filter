% Script to generate all needed plots.
% Plots include:
% - gamma comparison
% - Ensemble member comparison
% - Observation frequency comparison

% Given parameters
% initial and final time
t0 = 0;
tf = 10;
% number of time steps in solution
Nt = 10000;
% observation variance
sigma0 = sqrt(.01);

% Parameters of EnKPF_Standard
% number of observations
Nobs_standard = 200;
% ensemble members
Nens_standard = 2;
% gamma parameter in EnKPF
global H R gamma;
gamma_standard = 0.5;

% number of variables
Nvar = 3;

% initial condition
InitialCond = randn(Nvar,1);
% calculate time steps between observations
steps_standard = floor(Nt/Nobs_standard);
% setting up solution space
dt = (tf-t0)/Nt;
tSpace = linspace(t0,tf,Nt);

% generate the true solution
[T,XT] = SDESolver(dt,1,tf,InitialCond);
%XT = XT(2:end,:)';
%T = T(2:end);
TrueSolution = permute(XT,[3,1,2]);

% generate observations
XTObserved = TrueSolution(:,steps_standard:steps_standard:end);
tmp = eye(Nvar);
H = tmp(2:3,:);
%H = tmp;
R = sigma0*eye(size(H,1));
Observations = H*XTObserved + sigma0*randn(size(H,1),Nobs_standard);


% solution tracking memory allocation
EnKPF_standard = zeros(Nvar,Nens_standard,Nt);
EnKF_standard = zeros(Nvar,Nens_standard,Nt);
wPF_store = zeros(NP, numberOfSteps);
vPF_store = randn(varNum, NP, numberOfSteps);

% give each ensemble member a reasonable initial condition
EnKPF_standard(:,:,1) = bsxfun(@plus,InitialCond,0.3*randn(Nvar,Nens_standard));
EnKF_standard(:,:,1) = bsxfun(@plus,InitialCond,0.3*randn(Nvar,Nens_standard));

% EnKPF_standard
gamma = gamma_standard;
for ii=1:Nobs_standard
    disp(ii/Nobs_standard)
    
    % forecast ensemble
    [~,sol] = SDESolver(dt, Nens_standard, (steps_standard)*dt,...
        EnKPF_standard(:,:,(ii-1)*steps_standard+1)');
    
    EnKPF_standard(:,:,((ii-1)*steps_standard+1):(ii*steps_standard)) = ...
                permute(sol, [3,2,1]);
    
    % analysis update
    EnKPF_standard(:,:,ii*steps_standard) = ...
        EnKPF_update(EnKPF_standard(:,:,ii*steps_standard), Observations(:,ii), Nvar, Nens_standard);
    
    if (ii~=Nobs_standard)
        [~,sol] = SDESolver(dt, Nens_standard, dt, EnKPF_standard(:,:,ii*steps_standard)');
        EnKPF_standard(:,:,ii*steps_standard+1) = permute(sol, [3,2,1]);
    end
    
end

% EnKF_standard
for ii=1:Nobs_standard
    disp(ii/Nobs_standard)
    
    % forecaset ensemble
    [~,sol] = SDESolver(dt, Nens_standard, (steps_standard)*dt,...
        EnKF_standard(:,:,(ii-1)*steps_standard+1)');
    
    EnKF_standard(:,:,((ii-1)*steps_standard+1):(ii*steps_standard)) = ...
        permute(sol, [3,2,1]);
    
    % analysis update
    % compute the ensemble mean
    mu = (1/Nens_standard)*sum(EnKF_standard(:,:,ii*steps_standard), 2);
    
    % compute the ensemble covariance
    A = (bsxfun(@plus, EnKF_standard(:,:,ii*steps_standard), - mu))/(Nens_standard-1);
    
    % compute the Kalman gain matrix
    K = (A*(H*A)')/((H*A*(H*A)') + R);
    
    % apply the Kalman filter update on each ensemble member
    for jj=1:Nens_standard
        % Perturbation for y
        eps_y = normrnd(0,sigma0);
        EnKF_standard(:,jj,ii*steps_standard) = EnKF_standard(:,jj,ii*steps_standard)...
            + K*(Observations(:,ii) + eps_y - H*EnKF_standard(:,jj,ii*steps_standard));
    end
    
    if (ii~=Nobs_standard)
        [~,sol] = SDESolver(dt, Nens_standard, dt, EnKF_standard(:,:,ii*steps_standard)');
        EnKF_standard(:,:,ii*steps_standard+1) = permute(sol, [3,2,1]);
    end
    
    
end

% PF_standard
for ii = 1:Nobs_standard
    disp(ii/Nobs_standard)
    % PF
    % forecast particles
    [~,sol] = SDESolver(dt, Nens_standard, steps_standard*dt, vPF');
    vPF_store(:,:,((ii-1)*stepsBetweenObs+1):(ii*stepsBetweenObs)) = ...
                permute(sol(:,:,:), [3,2,1]);
    vPF = squeeze(sol(end,:,:))';
    
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

% plot the true solution and the ensemble mean at each time step
% only plot x, since that's the only variable that we care about
plot(tSpace,TrueSolution(1,:),...
    tSpace,permute(mean(EnKPF_standard(1,:,1:end),2),[3,2,1]),...
    tSpace,permute(mean(EnKF_standard(1,:,1:end),2),[3,2,1]));