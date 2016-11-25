% Implementation of the Ensemble Kalman Filter for the
% GBWB system. Generates the true solution using
% Euler-Maruyama, and tracks it with an ensemble kalman filter

ic = [0,0.2,0,0];
d = size(ic,2);
Nt = 1000;
NO = 10;
Nens = 2;
tfinal = 10;
obsVar = 0.3;

dt = (tfinal-ic(1))/Nt;
t = linspace(ic(1),tfinal,Nt);

trueSol = EulerMaruyama(@GBWB,ic,tfinal,Nt);
X = zeros(Nt,d,Nens);
Obs = zeros(NO,d);

for ii=1:Nens
    X(1,:,ii) = [ic(1),normrnd(0,0.3),normrnd(0,0.3),normrnd(0,0.3)];
end

for ii=2:Nt
    for jj=1:Nens
        X(ii,:,jj) = EulerMaruyama(@GBWB,X(ii-1,:,jj),X(ii,1,jj),1);
    end
    if (mod(ii,NO)==0)
        Obs(ii,:) = trueSol(ii,:) + normrnd(0,sqrt(obsVar));
        mu = (1/Nens)*sum(X(ii,:,:),3);
        A = reshape(X(ii,:,:),4,2) - mu;
        A = A(2:4,:);
        C = (A*A')/(Nens-1);
        K = C * (C + obsVar*eye(d))^-1;
        for jj=1:Nens
            X(ii,:,jj) = X(ii,:,jj) + K*(Obs(ii) - X(ii,:,jj));
        end
    end
end
plot(trueSol(:,1),trueSol(:,2),X(:,1,1),X(:,2,1));