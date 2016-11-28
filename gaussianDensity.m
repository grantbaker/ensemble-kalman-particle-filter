% This computes the gaussian density with mean mu and covariance C at x. k
% is the number of variables.
%
% This will be used to calculate the weights alpha.
%
function den = gaussianDensity(mu,C,x,k)
    den = ((2*pi)^k*det(C))*exp((-1/2)*(x-mu)'*(C\(x-mu)));
end