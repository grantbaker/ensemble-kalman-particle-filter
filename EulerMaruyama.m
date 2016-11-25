function [ sol ] = EulerMaruyama(stochDE, ic, tSpace)
%EULERMARUYAMA SDE solver
%   Self-explanatory

X = zeros(size(tSpace,2),size(ic,2));
X(1,:) = ic;

for ii=2:size(tSpace,2)
    dt = tSpace(ii)-tSpace(ii-1);
    X(ii,:) = X(ii-1,:) + stochDE(X(ii-1,:),1).*dt + stochDE(X(ii-1,:),2).*normrnd(0,sqrt(dt));
end

sol = X;


end

