% This function computes K(gamma*P_fore) according to formula (3) in Frei
% et al.
%
function K_gam_P = KalGain_gam_P(gamma, P_fore, H, R)
    K_gam_P = gamma*P_fore*H'/(gamma*H*P_fore*H' + R);
end