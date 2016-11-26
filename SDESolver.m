function [Time, Soln] = SDESolver(deltaT, NumSol, tFin)
    
    % Six tunable parameters
    a = 6;
    b = 3/4;
    c = -6;
    const = -3/4;
    gamma = 3;

    tInit = 0; % start time

    steps = floor((tFin - tInit)/deltaT); % Number of time steps
    
    Time = linspace(0, (steps-1)*deltaT, steps); %Initialize time vector
    
    Soln = zeros([steps NumSol 3]); % Initialize Vector of solutions
    % Initialize at one equilibrium
    Soln(1,:,:) = [1 0 0];
    
    % Only two variables are noisy.
    Noise = normrnd(0, sqrt(deltaT), [steps NumSol 2]); % Create noise increments
              
    % Note, First step of solution already set to zero

    % Solve SDE (Note all paths are being solved simultaneously)
    for k = 2:steps
        deltaX = (a*Soln(k-1,:,1) + b*(Soln(k-1,:,1)).^2 + ...
                  c*(Soln(k-1,:,1)).^3 + const + ...
                  gamma*Soln(k-1,:,2).*Soln(k-1,:,3))*deltaT;
        Soln(k,:,1) = Soln(k-1,:,1) + deltaX;

        deltaY = (-Soln(k-1,:,2) - gamma*Soln(k-1,:,1).*Soln(k-1,:,3))*deltaT + ...
                  sqrt(2)*Noise(k-1,:,1);
        Soln(k,:,2) = Soln(k-1,:,2) + deltaY;

        deltaZ = (-Soln(k-1,:,3))*deltaT + sqrt(2)*Noise(k-1,:,2);
        Soln(k,:,3) = Soln(k-1,:,3) + deltaZ;
    end
end
