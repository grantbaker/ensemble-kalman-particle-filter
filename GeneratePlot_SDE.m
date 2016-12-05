function GeneratePlot_SDE(TrueSolution, Nobs, Nens, gamma_, H_, t)
%GENERATEPLOT_SDE Filters using all three filtering methods, given the desired
%parameters. Generates a plot of the result.

global H R gamma Nvar Nt InitialCond sigma0 t0 tf dt;

gamma = gamma_;
H = H_;

% calculate time steps between observations
steps = floor(Nt/Nobs);

% generate observations
XTObserved = TrueSolution(:,steps:steps:end);

R = sigma0*eye(size(H,1));
Observations = H*XTObserved + sigma0*randn(size(H,1),Nobs);

% solution tracking memory allocation
EnKPF = zeros(Nvar,Nens,Nt);
EnKPFw = zeros(Nens,Nt);
EnKF = zeros(Nvar,Nens,Nt);
vPF = zeros(Nvar, Nens, Nt);
wPF = zeros(Nens,Nt);

% give each ensemble member a reasonable initial condition
EnKPF(:,:,1) = bsxfun(@plus,InitialCond,0.3*randn(Nvar,Nens));
EnKPFw(:,1:steps) = repmat(1/Nens,Nens,steps);
EnKF(:,:,1) = bsxfun(@plus,InitialCond,0.3*randn(Nvar,Nens));
vPF(:,:,1) = bsxfun(@plus,InitialCond,0.3*randn(Nvar,Nens));
wPF(:,1:steps) = repmat(1/Nens,Nens,steps);

% EnKPF_standard
for ii=1:Nobs
    disp(ii/Nobs)
    
    % forecast ensemble
    [~,sol] = SDESolver(dt, Nens, (steps)*dt,...
        EnKPF(:,:,(ii-1)*steps+1)');
    
    EnKPF(:,:,((ii-1)*steps+1):(ii*steps)) = ...
                permute(sol, [3,2,1]);
    
    % analysis update
    [update,alpha] = EnKPF_update(EnKPF(:,:,ii*steps), Observations(:,ii), Nvar, Nens);
    EnKPF(:,:,ii*steps) = update;
    EnKPFw(:,(ii*steps):((ii+1)*steps-1)) = repmat(alpha,1,steps);
    
    if (ii~=Nobs)
        [~,sol] = SDESolver(dt, Nens, dt, EnKPF(:,:,ii*steps)');
        EnKPF(:,:,ii*steps+1) = permute(sol, [3,2,1]);
    end
    
end

% EnKF_standard
for ii=1:Nobs
    disp(ii/Nobs)
    
    % forecaset ensemble
    [~,sol] = SDESolver(dt, Nens, (steps)*dt,...
        EnKF(:,:,(ii-1)*steps+1)');
    
    EnKF(:,:,((ii-1)*steps+1):(ii*steps)) = ...
        permute(sol, [3,2,1]);
    
    % analysis update
    % compute the ensemble mean
    mu = (1/Nens)*sum(EnKF(:,:,ii*steps), 2);
    
    % compute the ensemble covariance
    A = (bsxfun(@plus, EnKF(:,:,ii*steps), - mu))/sqrt(Nens-1);
    
    % compute the Kalman gain matrix
    K = (A*(H*A)')/((H*A*(H*A)') + R);
    
    % apply the Kalman filter update on each ensemble member
    for jj=1:Nens
        % Perturbation for y
        eps_y = normrnd(0,sigma0);
        EnKF(:,jj,ii*steps) = EnKF(:,jj,ii*steps)...
            + K*(Observations(:,ii) + eps_y - H*EnKF(:,jj,ii*steps));
    end
    
    if (ii~=Nobs)
        [~,sol] = SDESolver(dt, Nens, dt, EnKF(:,:,ii*steps)');
        EnKF(:,:,ii*steps+1) = permute(sol, [3,2,1]);
    end
    
    
end

% PF_standard
for ii = 1:Nobs
    disp(ii/Nobs)
    % PF
    % forecast particles
    [~,sol] = SDESolver(dt, Nens, steps*dt, vPF(:,:,(ii-1)*steps+1)');
    vPF(:,:,((ii-1)*steps+1):(ii*steps)) = ...
                permute(sol(:,:,:), [3,2,1]);
    
    for jj = 1:Nens
        wPF(jj,ii*steps) = exp(-.5*(norm(Observations(:,ii) - H*vPF(:,jj,(ii*steps-1)),2)/sigma0)^2);
    end
    
    wPF(:,(ii*steps):((ii+1)*steps-1)) = repmat(wPF(:,ii*steps)/sum(wPF(:,ii*steps)),1,steps);
    %wPF(:,((ii-1)*steps+1):(ii*steps)) = repmat(wPF(:,ii),steps,1)';

    % Resample
    NN = randsample(Nens,Nens,true,wPF(:,ii*steps));
    vPF(:,:,ii*steps) = vPF(:,NN,ii*steps);
    
    if (ii~=Nobs)
        [~,sol] = SDESolver(dt, Nens, dt, vPF(:,:,ii*steps)');
        vPF(:,:,ii*steps+1) = permute(sol, [3,2,1]);
    end
end

% plot the true solution and the ensemble mean at each time step
% only plot x, since that's the only variable that we care about

tSpace = linspace(t0,tf,Nt);

mean_EnKPF = permute(squeeze(sum(bsxfun(@times, permute(EnKPF, [2,3,1]), EnKPFw(:,1:Nt)),1)),[2,1]);
mean_PF = permute(squeeze(sum(bsxfun(@times, permute(vPF, [2,3,1]), wPF(:,1:Nt)),1)),[2,1]);

figure;
set(0,'defaultaxesfontname','courier');
set(0,'defaulttextinterpreter','latex');
set(0, 'defaultLegendInterpreter','latex')
p=plot(tSpace,mean_EnKPF(1,:,:),...
    tSpace,permute(mean(EnKF(1,:,1:end),2),[3,2,1]),...
    tSpace,mean_PF(1,:,:),...
    tSpace,TrueSolution(1,:));
p(4).LineWidth = 2;
legend('EnKPF','EnKF','PF','True')
%set(a,'TickLabelInterpreter', 'latex');
title(t)

RMS_EnKPF = sqrt(mean((mean_EnKPF-TrueSolution).^2));
RMS_EnKF = sqrt(mean((squeeze(mean(EnKF,2))-TrueSolution).^2));
RMS_PF = sqrt(mean((mean_PF-TrueSolution).^2));

figure;
set(0,'defaultaxesfontname','courier');
set(0,'defaulttextinterpreter','latex');
set(0, 'defaultLegendInterpreter','latex')
p=plot(tSpace,RMS_EnKPF,...
    tSpace,RMS_EnKF,...
    tSpace,RMS_PF);
legend('EnKPF','EnKF','PF')
title(t)


end