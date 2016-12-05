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
sigma0 = sqrt(.01);

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

GeneratePlot_SDE(TrueSolution, 200, 30, 0.5, H1, 'Standard')
GeneratePlot_SDE(TrueSolution, 200, 30, 0.5, H1, 'Standard')
GeneratePlot_SDE(TrueSolution, 200, 30, 0.5, H2, 'Observe $x$')
GeneratePlot_SDE(TrueSolution, 200, 30, 0.5, H2, 'Observe $x$')

GeneratePlot_SDE(TrueSolution, 200, 30, 0.05, H1, '$\gamma = 0.05$')
GeneratePlot_SDE(TrueSolution, 200, 30, 0.05, H1, '$\gamma = 0.05$')
GeneratePlot_SDE(TrueSolution, 200, 30, 0.3, H1, '$\gamma = 0.3$')
GeneratePlot_SDE(TrueSolution, 200, 30, 0.3, H1, '$\gamma = 0.3$')
GeneratePlot_SDE(TrueSolution, 200, 30, 0.7, H1, '$\gamma = 0.7$')
GeneratePlot_SDE(TrueSolution, 200, 30, 0.7, H1, '$\gamma = 0.7$')
GeneratePlot_SDE(TrueSolution, 200, 30, 0.95, H1, '$\gamma = 0.95$')
GeneratePlot_SDE(TrueSolution, 200, 30, 0.95, H1, '$\gamma = 0.95$')

GeneratePlot_SDE(TrueSolution, 20, 30, 0.5, H1, '20 Observations')
GeneratePlot_SDE(TrueSolution, 20, 30, 0.5, H1, '20 Observations')
GeneratePlot_SDE(TrueSolution, 2000, 30, 0.5, H1, '2000 Observations')
GeneratePlot_SDE(TrueSolution, 2000, 30, 0.5, H1, '2000 Observations')

GeneratePlot_SDE(TrueSolution, 200, 10, 0.5, H1, '10 Ensemble Members')
GeneratePlot_SDE(TrueSolution, 200, 10, 0.5, H1, '10 Ensemble Members')
GeneratePlot_SDE(TrueSolution, 200, 100, 0.5, H1, '100 Ensemble Members')
GeneratePlot_SDE(TrueSolution, 200, 100, 0.5, H1, '100 Ensemble Members')
GeneratePlot_SDE(TrueSolution, 200, 500, 0.5, H1, '500 Ensemble Members')
GeneratePlot_SDE(TrueSolution, 200, 500, 0.5, H1, '500 Ensemble Members')