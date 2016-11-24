function [ dX ] = GBWB(Xt,stoch)
%GBWB A three dimensional system of SDEs
%   more stuff

a=-6;
b=-3/4;
c=6;
k=3/4;
gamma=3;

if (stoch==1)
    dX=[1,-(a*Xt(2) + b*Xt(2)^2 + c*Xt(2)^3 + k) + gamma*Xt(3)*Xt(4),-Xt(3) - gamma*Xt(2)*Xt(4), -Xt(4)];
else
    dX=[0,0,sqrt(2),sqrt(2)];
end

end

