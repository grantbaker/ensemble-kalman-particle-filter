% Implementation of the Ensemble Kalman Filter for the
% GBWB system. Generates the true solution using
% Euler-Maruyama, and tracks it with an ensemble kalman filter

trueSol = EulerMaruyama(@GBWB,[0,0.2,0,0],10,100000);

