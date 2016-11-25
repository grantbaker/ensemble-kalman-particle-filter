% Implementation of the Ensemble Kalman Filter for the
% GBWB system. Generates the true solution using
% Euler-Maruyama, and tracks it with an ensemble kalman filter

ic = [0.2,0,0];
d = size(ic,2);
t0 = 0;
tf = 10;
Nt = 10000;
NO = 1;
Nens = 100;
obsVar = 0.0001;

dt = (tf-t0)/Nt;
tSpace = linspace(t0,tf,Nt);

trueSol = EulerMaruyama(@GBWB,ic,tSpace);
X = zeros(Nt,d,Nens);
Obs = zeros(NO,d);

for ii=1:Nens
    X(1,:,ii) = [normrnd(0,0.3),normrnd(0,0.3),normrnd(0,0.3)];
end

for ii=2:Nt
    for jj=1:Nens
        tmp = EulerMaruyama(@GBWB,X(ii-1,:,jj),[tSpace(ii-1),tSpace(ii)]);
        X(ii,:,jj) = tmp(2,:);
    end
    if (mod(ii,NO)==0)
        Obs(ii,:) = trueSol(ii,:) + normrnd(0,sqrt(obsVar));
        mu = (1/Nens)*sum(X(ii,:,:),3);
        A = reshape(X(ii,:,:),3,Nens) - mu';
        C = (A*A')/(Nens-1);
        K = C * (C + obsVar*eye(d))^-1;
        for jj=1:Nens
            X(ii,:,jj) = (X(ii,:,jj)' + K*(Obs(ii)' - X(ii,:,jj)'))';
        end
    end
end
plot(tSpace,trueSol(:,1),tSpace,mean(X(:,1,:),3));