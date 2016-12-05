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
Nvar = 3;
% observation matrix
tmp = eye(Nvar);
H1 = tmp(2:3,:);
H2 = tmp;

% initial condition
InitialCond = randn(Nvar,1);
% setting up solution space
dt = (tf-t0)/Nt;

% generate the true solution
[T,XT] = SDESolver(dt,1,tf,InitialCond);
%XT = XT(2:end,:)';
%T = T(2:end);
TrueSolution = permute(XT,[3,1,2]);

GeneratePlot(TrueSolution, 500, 200, 0.5, H1)
GeneratePlot(TrueSolution, 500, 200, 0.5, H2)