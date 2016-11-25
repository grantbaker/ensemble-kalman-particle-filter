function [ dX ] = brownian(~, stoch )
%BROWNIAN brownian motion

if (stoch==1)
    dX=0;
else
    dX=1;

end

