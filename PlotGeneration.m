% Script to generate all needed plots.
% Plots include:
% - gamma comparison
% - Ensemble member comparison
% - Observation frequency comparison

% Given Constan
% initial and final time
t0 = 0;
tf = 10;
% number of time steps in solution
Nt = 10000;
% number of observations
Nobs = 200;
% ensemble members
Nens = 20;
% observation variance
sigma0 = sqrt(.2);
% gamma parameter in EnKPF
global H R gamma;
gamma = 0.5;

% number of variables
Nvar = 3;

% initial condition
InitialCond = randn(Nvar,1);
% calculate time steps between observations
stepsBetweenObs = floor(Nt/Nobs);
% setting up dimension of system and solution space
dt = (tf-t0)/Nt;
tSpace = linspace(t0,tf,Nt);

% generate the true solution
[T,XT] = SDESolver(dt,1,tf,InitialCond);
%XT = XT(2:end,:)';
%T = T(2:end);
TrueSolution = permute(XT,[3,1,2]);

% generate observations
XTObserved = TrueSolution(:,stepsBetweenObs:stepsBetweenObs:end);
tmp = eye(Nvar);
%H = tmp(2:3,:);
H = tmp;
R = sigma0*eye(size(H,1));
Observations = H*XTObserved + sigma0*randn(size(H,1),Nobs);


% solution tracking memory allocation
EnKPF_standard = zeros(Nvar,Nens,Nt);

% give each ensemble member a reasonable initial condition
EnKPF_standard(:,:,1) = bsxfun(@plus,InitialCond,0.3*randn(Nvar,Nens));

% begin propogating ensemble
for ii=1:Nobs
    disp(ii/Nobs)
    
    % EnKPF_standard
    % forecast ensemble
    [~,sol] = SDESolver(dt, Nens, (stepsBetweenObs)*dt,...
        EnKPF_standard(:,:,(ii-1)*stepsBetweenObs+1)');
    
    EnKPF_standard(:,:,((ii-1)*stepsBetweenObs+1):(ii*stepsBetweenObs)) = ...
                permute(sol, [3,2,1]);
    
    % analysis update
    EnKPF_standard(:,:,ii*stepsBetweenObs) = ...
        EnKPF_update(EnKPF_standard(:,:,ii*stepsBetweenObs), Observations(:,ii), Nvar, Nens);
    
    if (ii~=Nobs)
        [~,sol] = SDESolver(dt, Nens, dt, EnKPF_standard(:,:,ii*stepsBetweenObs)');
        EnKPF_standard(:,:,ii*stepsBetweenObs+1) = permute(sol, [3,2,1]);
    end
    
end

% plot the true solution and the ensemble mean at each time step
% only plot x, since that's the only variable that we care about
plot(tSpace,TrueSolution(1,:),tSpace,permute(mean(EnKPF_standard(1,:,1:end),2),[3,2,1]));