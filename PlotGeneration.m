% Script to generate all needed plots.
% Plots include:
% - gamma comparison
% - Ensemble member comparison
% - Observation frequency comparison


% Parameters
% initial condition
ic = [0.2,0,0];
% initial and final time
t0 = 0;
tf = 10;
% number of time steps in solution
Nt = 10000;
% number of time steps between assimilations
NO = 2000;
% ensemble members
Nens = 20;
% observation variance
obsVar = 0.3;

% setting up dimension of system and solution space
d = size(ic,2);
dt = (tf-t0)/Nt;
tSpace = linspace(t0,tf,Nt);

% generate the true solution
trueSol = EulerMaruyama(@GBWB,ic,tSpace);
% allocate memory
X = zeros(Nt,d,Nens);
Obs = zeros(NO,d);

% give each ensemnble member a reasonable initial condition
for ii=1:Nens
    X(1,:,ii) = [normrnd(0,0.3),normrnd(0,0.3),normrnd(0,0.3)];
end

% begin propogating ensemble
for ii=2:Nt
    %disp(ii/Nt);
    % solve each ensemble member to the next time
    for jj=1:Nens
        tmp = EulerMaruyama(@GBWB,X(ii-1,:,jj),[tSpace(ii-1),tSpace(ii)]);
        X(ii,:,jj) = tmp(2,:);
    end
    
    % if adequate number of steps has passed, assimilate
    if (mod(ii,NO)==0)
        disp(ii/Nt);
        
        % generate observation
	    Obs(ii,:) = trueSol(ii,:) + normrnd(0,sqrt(obsVar));
        
        % compute the ensemble mean
        mu = (1/Nens)*sum(X(ii,:,:),d);
        
        % compute the ensemble covariance
        %A = reshape(X(ii,:,:),d,Nens)- mu';
        A = reshape(X(ii,:,:),d,Nens);
        for jj=1:Nens
            A(:,jj) = A(:,jj)-mu';
        end
        C = (A*A')/(Nens-1);
        
        % compute the Kalman gain matrix
        K = C * (C + obsVar*eye(d))^-1;
        
        % apply the Kalman filter update on each ensemble member
        for jj=1:Nens
            X(ii,:,jj) = (X(ii,:,jj)' + K*(Obs(ii)' - X(ii,:,jj)'))';
        end
    end
end

% plot the true solution and the ensemble mean at each time step
% only plot x, since that's the only variable that we care about
plot(tSpace,trueSol(:,1),tSpace,mean(X(:,1,:),3));