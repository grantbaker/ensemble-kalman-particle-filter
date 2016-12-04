function [ dX ] = GBWB(Xt,stoch)
%GBWB A three dimensional system of SDEs
%   more stuff

% proposed parameters
%a=-6;
%b=-3/4;
%c=6;
%k=3/4;
%gamma=1;

% true parameters
a=6;
b=3/4;
c=-6;
k=-3/4;
gamma=3;

if (stoch==1)
    dX=[a*Xt(1) + b*Xt(1)^2 + c*Xt(1)^3 + k + gamma*Xt(2)*Xt(3),-Xt(2) - gamma*Xt(1)*Xt(3), -Xt(3)];
else
    dX=[0,sqrt(2),sqrt(2)];
end

end

