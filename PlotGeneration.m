% Script to generate all needed plots.
% Plots include:
% - gamma comparison
% - Ensemble member comparison
% - Observation frequency comparison

% Parameters
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
obsVar = 0.2;

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
[T,trueSolution] = SDESolver(dt,1,tf,InitialCond);

% generate observations
XTObserved = XT(:,1:stepsBetweenObs:end);
tmp = eye(varNum);
H = tmp(2:3,:);
%H = tmp;
R = obsVar*eye(size(H,1));
Observations = H*XTObserved + obsVar*randn(size(H,1),200);


% give each ensemble member a reasonable initial condition
EnKPF = bsxfun(@plus,InitialCond,0.3*randn(Nvar,Nens));

% solution tracking memory allocation
EnKPF_solution = zeros(Nt,Nvar,Nens);

% begin propogating ensemble
for ii=1:Nobs
    %disp(ii/Nt);
    % solve each ensemble member to the next time
    
    % if adequate number of steps has passed, assimilate
    if (mod(ii,Nobs)==0)
        disp(ii/Nt);
        
        % generate observation
	    Observations(ii,:) = trueSol(ii,:) + normrnd(0,sqrt(obsVar));
        
        % compute the ensemble mean
        mu = (1/Nens)*sum(X(ii,:,:),Nvar);
        
        % compute the ensemble covariance
        %A = reshape(X(ii,:,:),d,Nens)- mu';
        A = reshape(X(ii,:,:),Nvar,Nens);
        for jj=1:Nens
            A(:,jj) = A(:,jj)-mu';
        end
        C = (A*A')/(Nens-1);
        
        % compute the Kalman gain matrix
        K = C * (C + obsVar*eye(Nvar))^-1;
        
        % apply the Kalman filter update on each ensemble member
        for jj=1:Nens
            X(ii,:,jj) = (X(ii,:,jj)' + K*(Observations(ii)' - X(ii,:,jj)'))';
        end
    end
end

% plot the true solution and the ensemble mean at each time step
% only plot x, since that's the only variable that we care about
plot(tSpace,trueSol(:,1),tSpace,mean(X(:,1,:),3));