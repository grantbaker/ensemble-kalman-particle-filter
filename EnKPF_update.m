% This function runs one update step of the EnKPF described in Frei et al.
% 2012.
%
% Inputs: fore is the ensemble forecast (this is varNum by ensNum), 
%         obs is the vector of observations (this is d by 1), 
%         varNum is the dimension of the system, 
%         ensNum is the number of ensemble members,
%
% Global: gamma is a parameter in [0,1].
%         H is Observation matrix (d by varNum),
%         R is Observation Error covariance (d by d).
%
% Output: the ensemble update update (this is varNum x ensNum)
%
function update = EnKPF_update(fore, obs, varNum, ensNum)

    global gamma H R;
                            
    % Initiate matrices
    update = zeros(varNum, ensNum);
    nu = zeros(varNum, ensNum);
    alpha = zeros(ensNum, 1);

    % 1. Compute the estimated forecast covariance P_fore.
    % - Use built in unbiased covariance to estimate.

    P_fore = cov(fore',0);

    % 2. Choose gamma and compute K(gamma*P_fore) and nu.
    % - Use formulas (3) and (5) in Frei et al.

    % gamma is given in the function arguments.
    K_gam_P = KalGain_gam_P(gamma, P_fore, H, R);
    Q_gam_P = compute_QgamP(gamma, K_gam_P, R);

    % 3. Compute Q(gamma, P_fore) and wieghts, alpha.
    % - Use formulas (6) and (8) in Frei et al.

    for ii = 1:ensNum
        x_prior = fore(:, ii);
        % Compute the nu values
        nu(:, ii) = compute_nu(x_prior, K_gam_P, obs, H);
        % Compute the weights. Will need to normalize out of the loop.
        alpha(ii) = gaussianDensity(H*nu(:,ii), ...
                                    H*Q_gam_P*(H') + (1/(1-gamma))*R, ...
                                    obs, varNum);
    end
    alpha = alpha/(sum(alpha));

    % 4. Choose indices I(j) by sampling from the weights alpha with the
    % balanced sampling scheme in Kuensch (2005).
    % - Use built in randsample from matlab.

    Indicies = randsample(ensNum, ensNum, true, alpha);

    % 5. Generate err_1 ~ N(0,R) and set 
    % x_gam(j) = nu(I(j)) + K(gamma*P_fore)*(err_1(j)/sqrt(gamma))

    % 6. Compute K((1-gamma)Q(gamma,P_fore)), generate err_2 ~ N(0,R) and set
    % x_update(j) = x_gam(j) + K((1-gamma)Q(gamma,P_fore))(obs + ...
    % (eps_2(j)/sqrt(1-gamma)) - Hx_gam(j)).

    KalGain2 = KalGain_gam_P(1-gamma, Q_gam_P, H, R);

    for ii = 1:ensNum
        
        err_1 = mvnrnd(zeros(size(R,1),1), R)';
        err_2 = mvnrnd(zeros(size(R,1),1), R)';

        x_gam = nu(:,Indicies(ii)) + K_gam_P*(err_1/sqrt(gamma));
        update(:,ii) = x_gam + KalGain2*(obs + (err_2/sqrt(1-gamma)) - ...
                        H*x_gam);
    end

end