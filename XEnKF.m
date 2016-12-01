% Implementation of the Mixture Ensemble Kalman Filter for the
% GBWB system. Generates the true solution using
% Euler-Maruyama, and tracks it with an ensemble kalman filter
% as well as a mixture ensemble kalman filter

% Parameters
% initial condition
ic = [0.2,0,0];
% initial and final time
t0 = 0;
tf = 1;
% number of time steps in solution
Nt = 1000;
% number of time steps between assimilations
NO = 10;
% number of fuzzy clusters for XEnKF
Nl = 10;
% ensemble members
Nens = 100;
% observation variance
obsVar = 0.2;

% setting up dimension of system and solution space
d = size(ic,2);
dt = (tf-t0)/Nt;
tSpace = linspace(t0,tf,Nt);

% fcm clustering method options
fcmoptions = [2.0, 200, 1e-5, 1];

% generate the true solution
trueSol = EulerMaruyama(@GBWB,ic,tSpace);
clear ic;

% allocate memory
X = zeros(Nt,d,Nens);
Xl = zeros(Nt,d,Nens);
Obs = zeros(NO,d);

% give each ensemnble member a reasonable initial condition
for ii=1:Nens
    ictmp = [normrnd(0,0.3),normrnd(0,0.3),normrnd(0,0.3)];
    % the ensemble members start out with the same ICs
    X(1,:,ii) = ictmp;
    Xl(1,:,ii) = ictmp;
end
clear ictmp;

% begin propogating ensemble and updating
for ii=2:Nt
    %disp(ii/Nt);
    % solve each ensemble member to the next time
    for jj=1:Nens
        % propogation for EnKF ensemble
        tmp = EulerMaruyama(@GBWB,X(ii-1,:,jj),[tSpace(ii-1),tSpace(ii)]);
        X(ii,:,jj) = tmp(2,:);
        % propogation for XEnKF ensemble
        tmp = EulerMaruyama(@GBWB,Xl(ii-1,:,jj),[tSpace(ii-1),tSpace(ii)]);
        Xl(ii,:,jj) = tmp(2,:);
    end
    
    % if adequate number of steps has passed, assimilate
    if (mod(ii,NO)==0)
        %disp(ii/Nt);
        
        % generate observation
	    Obs(ii,:) = trueSol(ii,:) + normrnd(0,sqrt(obsVar));
        
        % ---- ENKF UPDATE ----
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
            X(ii,:,jj) = (X(ii,:,jj)' + K*(Obs(ii,:)' - X(ii,:,jj)'))';
        end
        
        % ---- XENKF UPDATE ----
        
        % applying the fuzzy c-means clustering method
        [centers,psi] = fcm(reshape(Xl(ii,:,:),d,Nens)',Nl,fcmoptions);
        
        % compute alpha for each cluster
        alpha = (1/Nens)*sum(psi,2);
        
        % compute tau for each cluster
        tau = zeros(size(psi));
        for jj=1:Nl
            tau(jj,:) = psi(jj,:)/(Nens*alpha(jj));
        end
        
        % compute the mean and covariance
        muX = zeros(d,Nl);
        for jj=1:Nl
            for k=1:Nens
                muX(:,jj) = muX(:,jj) + tau(jj,k)*Xl(ii,:,k)';
            end
        end
        P = zeros(d,d,Nl);
        for jj=1:Nl
            for k=1:Nens
                P(:,:,jj) = P(:,:,jj) + tau(jj,k)*(Xl(ii,:,k)'-muX(:,jj))*(Xl(ii,:,k)'-muX(:,jj))';
            end
        end
        
        % compute the updated weights alphaU and normalize
        alphaU = zeros(size(alpha));
        for jj=1:Nl
            alphaU(jj) = alpha(jj)*mvnpdf(Obs(ii,:)',muX(:,jj),(P(:,:,jj)+(obsVar)*eye(d)));
        end
        alphaU = alphaU/sum(alphaU);
        
        % generate index pairs
        indices = zeros(Nens,2);
        indices(:,1) = randsample(1:Nl,Nens,true,alphaU);
        for jj=1:Nens
            indices(jj,2) = randsample(1:Nens,1,true,tau(indices(jj,1),:)');
        end
        
        % update ensemble
        for jj=1:Nens
            % compute kalman gain
            Kx = P(:,:,indices(jj,1))*(P(:,:,indices(jj,1)) + obsVar*eye(d))^-1;
            Xl(ii,:,jj) = (Xl(ii,:,indices(jj,2))' + Kx*(Obs(ii,:)' + normrnd(0,sqrt(obsVar)) - Xl(ii,indices(jj,2))'))';
        end
        
    end
end

% plot the true solution and the ensemble mean at each time step
% only plot x, since that's the only variable that we care about
plot(tSpace,trueSol(:,1),tSpace,mean(X(:,1,:),3),tSpace,mean(Xl(:,1,:),3));
