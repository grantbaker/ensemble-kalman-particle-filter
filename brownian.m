function [ dX ] = brownian(Xt, stoch )
%BROWNIAN brownian motion

if (stoch==1)
    dX=[1,0];
else
    dX=[0,1];

end

