function [ sol ] = EulerMaruyama(stochDE, X0, t1, N)
%EULERMARUYAMA SDE solver
%   Self-explanator

X = zeros(N,size(X0));
X(1) = X0;

dt = (t1-X0(1))/N;

for ii=2:N
    X(ii) = X(ii-1) + stochDE(X(ii-1),1)*dt + stochDE(X(ii-1),2)*normrnd(0,dt);
end

sol = X;


end

