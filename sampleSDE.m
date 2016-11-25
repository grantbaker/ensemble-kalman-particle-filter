function [ dX ] = sampleSDE(Xt, stoch)
%SAMPLESDE a sample SDE to solve.
%   The Ornstein-Uhlenbeck process

% parameters
theta = 0.1;
mu = 1.5;
sigma = 6;

if (stoch==1)
    dX = theta*(mu-Xt);
else
    dX = sigma;
end


end

