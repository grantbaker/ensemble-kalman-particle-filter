% This function computes Q(gamma, P_fore) according to formula (5) in Frei
% et al.
%
function Q_gam_P = compute_QgamP(gamma, KalGain, R)
    Q_gam_P = (1/gamma)*KalGain*R*(KalGain');
end