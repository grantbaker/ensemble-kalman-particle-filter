function [ sol ] = EulerMaruyama(stochDE, ic, tSpace)
%EULERMARUYAMA SDE solver using the Euler-Maruyama method
%   Solves the system of stochastic differential equations given by 
% StochDE, starting at the initial condition ic and using each value in 
% tSpace as a point. Outputs the solution sol at each time in tSpace.
% Typically, tSpace is a linspace from some initial time to a final time.

% allocate memory
X = zeros(size(tSpace,2),size(ic,2));
X(1,:) = ic;

% solve SDE at each point in tSpace using Euler-Maruyama
for ii=2:size(tSpace,2)
    dt = tSpace(ii)-tSpace(ii-1);
    X(ii,:) = X(ii-1,:) + stochDE(X(ii-1,:),1).*dt + stochDE(X(ii-1,:),2).*normrnd(0,sqrt(dt));
end

sol = X;


end

