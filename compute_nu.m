% This function computes nu for a single ensemble member according to
% formula (5) from Frei et al.
%
function nu = compute_nu(x_prior, KalGain, obs, H)
    nu = x_prior + KalGain*(obs - H*x_prior);
end