% This implements one step of the XEnKF update
%
% Inputs: fore is the ensemble forecast (this is varNum by ensNum), 
%         obs is the vector of observations (this is d by 1), 
%         prev_mixingProb is the previous vector of mixing probabilities
%               (this is numCenters by 1)
%         varNum is the dimension of the system, 
%         numCenters is the number of centers for the gaussian mixture, 
%         N is the number of nearest neighbors to each center, and 
%         ensNum is the number of ensemble members
%
% Note: varNum >= numCenters >= N
%
% Global: Observation matrix H is (d by varNum)
%         Observation Error R covariance is (d by d)
%
% Output: the ensemble update update (this is varNum x ensNum)
function [update, mixingProb] = XEnKF_update(fore, obs, prev_mixingProb, ...
                                varNum, numCenters, N, ensNum)
    
    global H R;
    
    % Because it does not matter which variables we pick to be the centers,
    % we pick 1 to numCenters.
    update = zeros(varNum,ensNum);
                            
    % Matrix of N nearest neighbors for each center.
        % The rows correspond to the ordered list of centers for the
        % gaussian mixture, and the columns correspond to the N nearest
        % neighbors.
    nearestN = zeros(numCenters,N);
    
    for l = 1:numCenters
        % Find the nearest N neighbors of each center.
        indices_for_lth_center = knnsearch((fore(:, [1:(l-1) (l+1):end]))', ...
                                            (fore(:,l))', 'K', N);
                     
        % Note, the above search finds the index of the item in the list
        % fore(:, [1:(l-1) (l+1):end]) which is one shorted than the entire
        % list fore. So, if we are beyond l, the index is one off and must
        % be shifted.
        for k = 1:N
            if (indices_for_lth_center(k) >= l)
                indices_for_lth_center(k) = indices_for_lth_center(k) + 1;
            end
        end
                                        
        nearestN(l,:) = indices_for_lth_center;
    end
    
    % Find values to calculate
    wVec = zeros(numCenters,1);
    for k = 1:numCenters
        % Obtain the indices of the nearest vectors to the current center
        nearestVec = squeeze(nearestN(k,:));
        % Obtain the varNum by N matrix of data
        data = fore(:,nearestVec);
        % Input of cov assumes columns are variables and rows are
        % observations. Thus, we transpose data.
        P_centers = cov(data',0);
        % Now calculate the mean for each center
        mu_centers = mean(data, 2);
        
        % Calculate mixing probabilities
        wVec(k) = ((det(H*P_centers*(H') + R))^(-1/2))* ...
                        exp(-(1/2)*((obs - H*mu_centers)')* ...
                        inv(H*P_centers*(H') + R)* ...
                        (obs - H*mu_centers));
    end
    
    % Finally, calculate the mixing probabilities.
    mixingProb = prev_mixingProb.*wVec;
    mixingProb = mixingProb/sum(mixingProb);
    
    for j = 1:ensNum
        % Index of xI
        Ind = randsample(numCenters,1,true,mixingProb);
        % Random Index of point near xI;
        Ind2 = unidrnd(N);
        
        % Find the correct nearest neighbor of xI
        xI = fore(:, Ind);
        nearestVec = squeeze(nearestN(Ind,:));
        data_near_xI = fore(:,nearestVec);
        
        indexOfxstar = nearestN(Ind,Ind2);
        xstar = fore(:, indexOfxstar);
        
        % Need to invert data for cov to work right
        covI = cov(data_near_xI',0);
        
        kalGain = covI*(H')*inv(H*covI*(H') + R);
        
        err = mvnrnd(zeros(size(R,1),1),R)';
        
        update(:,j) = xstar + kalGain*(obs + err - H*xstar);
    end
end
