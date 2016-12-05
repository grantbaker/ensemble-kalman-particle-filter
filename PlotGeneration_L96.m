% Script to generate all needed plots.
% Plots include:
% - gamma comparison
% - Ensemble member comparison
% - Observation frequency comparison

global H R gamma Nvar Nt InitialCond sigma0 tf t0 dt;

% Given parameters
% initial and final time
t0 = 0;
tf = 10;
% number of time steps in solution
Nt = 10000;
% observation variance
sigma0 = sqrt(.2);

% number of variables
Nvar = 40;
% observation matrix
tmp = eye(Nvar);
H1 = tmp(2:2:end,:);

% initial condition
InitialCond = randn(Nvar,1);
% setting up solution space
dt = (tf-t0)/Nt;

% generate the true solution
[T,XT] = ode45(@RHS_L96,linspace(t0,tf,Nt),randn(Nvar,1));
%XT = XT(2:end,:)';
%T = T(2:end);
TrueSolution = permute(XT,[2,1]);

GeneratePlot_L96(TrueSolution, 500, 100, 0.5, H1, 'Standard')