function GeneratePlot_L96(TrueSolution, Nobs, Nens, gamma_, H_, t)
%GENERATEPLOT_L96 Filters using all three filtering methods, given the desired
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
wPF = zeros(Nens,Nt);
vPF = zeros(Nvar, Nens, Nt);

% give each ensemble member a reasonable initial condition
EnKPF(:,:,1) = bsxfun(@plus,InitialCond,0.3*randn(Nvar,Nens));
EnKPFw(:,1:steps) = repmat(1/Nens,Nens,steps);
EnKF(:,:,1) = bsxfun(@plus,InitialCond,0.3*randn(Nvar,Nens));
%vPF(:,:,1) = bsxfun(@plus,InitialCond,0.3*randn(Nvar,Nens));
%wPF(:,1:steps) = repmat(1/Nens,Nens,steps);

% EnKPF_standard
for ii=1:Nobs
    disp(ii/Nobs)
    
    % forecast ensemble
    for k=1:Nens
        [~,sol] = ode45(@RHS_L96, 0:dt:((steps-1)*dt), EnKPF(:,k,(ii-1)*steps+1));
    
        EnKPF(:,k,((ii-1)*steps+1):(ii*steps)) = permute(sol, [2,3,1]);
    end
    
    % analysis update
    [update,alpha] = EnKPF_update(EnKPF(:,:,ii*steps), Observations(:,ii), Nvar, Nens);
    EnKPF(:,:,ii*steps) = update;
    EnKPFw(:,(ii*steps):((ii+1)*steps-1)) = repmat(alpha,1,steps);
    
    if (ii~=Nobs)
        for k=1:Nens
            [~,sol] = ode45(@RHS_L96, [0 dt], EnKPF(:,k,ii*steps));
            EnKPF(:,k,ii*steps+1) = permute(sol(end,:), [2,3,1]);
        end
    end
    
end

% EnKF_standard
for ii=1:Nobs
    disp(ii/Nobs)
    
    % forecaset ensemble
    for k=1:Nens
        [~,sol] = ode45(@RHS_L96, 0:dt:((steps-1)*dt), EnKF(:,k,(ii-1)*steps+1));
    
        EnKF(:,k,((ii-1)*steps+1):(ii*steps)) = permute(sol, [2,3,1]);
    end
    
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
        for k=1:Nens
            [~,sol] = ode45(@RHS_L96, [0 dt], EnKF(:,k,ii*steps));
            EnKF(:,k,ii*steps+1) = permute(sol(end,:), [2,3,1]);
        end
    end
    
    
end

% PF_standard
% for ii = 1:Nobs
%     disp(ii/Nobs)
%     % PF
%     % forecast particles
%     for k=1:Nens
%         [~,sol] = ode45(@RHS_L96, 0:dt:((steps-1)*dt), vPF(:,k,(ii-1)*steps+1));
%     
%         vPF(:,k,((ii-1)*steps+1):(ii*steps)) = permute(sol, [2,3,1]);
%     end
%     
%     for jj = 1:Nens
%         wPF(jj) = exp(-.5*(norm(Observations(:,ii) - H*vPF(:,jj,ii*steps),2)/sigma0)^2);
%     end
%     
%     wPF = wPF/sum(wPF);
%     %wPF(:,((ii-1)*steps+1):(ii*steps)) = repmat(wPF(:,ii),steps,1)';
% 
%     % Resample
%     NN = randsample(Nens,Nens,true,wPF);
%     vPF(:,:,ii*steps) = vPF(:,NN,ii*steps);
%     
%     if (ii~=Nobs)
%         for k=1:Nens
%             [~,sol] = ode45(@RHS_L96, [0 dt], vPF(:,k,ii*steps));
%             vPF(:,k,ii*steps+1) = permute(sol(end,:), [2,3,1]);
%         end
%     end
% end

% plot the true solution and the ensemble mean at each time step
% only plot x, since that's the only variable that we care about
mean_EnKPF = permute(squeeze(sum(bsxfun(@times, permute(EnKPF, [2,3,1]), EnKPFw(:,1:Nt)),1)),[2,1]);

tSpace = linspace(t0,tf,Nt);

% figure;
% set(0,'defaultaxesfontname','courier');
% set(0,'defaulttextinterpreter','latex');
% set(0, 'defaultLegendInterpreter','latex')
% p=plot(tSpace,permute(mean(EnKPF(1,:,1:end),2),[3,2,1]),...
%     tSpace,permute(mean(EnKF(1,:,1:end),2),[3,2,1]),...
%     tSpace,permute(mean(vPF(1,:,1:end),2),[3,2,1]),...
%     tSpace,TrueSolution(1,:));
% p(4).LineWidth = 2;
% legend('EnKPF','EnKF','PF','True')
% %set(a,'TickLabelInterpreter', 'latex');
% title(t)

RMS_EnKPF = sqrt(mean((mean_EnKPF-TrueSolution).^2));
RMS_EnKF = sqrt(mean((squeeze(mean(EnKF,2))-TrueSolution).^2));

RMS_EnKPF_mean = mean(RMS_EnKPF);
RMS_EnKF_mean = mean(RMS_EnKF);

%RMS_PF = sqrt(mean((squeeze(mean(vPF,2))-TrueSolution).^2));
%tSpace,RMS_PF
figure;
set(0,'defaultaxesfontname','courier');
set(0,'defaulttextinterpreter','latex');
set(0, 'defaultLegendInterpreter','latex')
p=plot(tSpace,RMS_EnKPF,...
    tSpace,RMS_EnKF);
leg = legend(['EnKPF; RMS = ',num2str(RMS_EnKPF_mean)],['EnKF; RMS = ',num2str(RMS_EnKF_mean)]);
legtxt = findobj(leg, 'type', 'text');
set(legtxt,'FontSize',18)
title(t,'FontSize',18)


end

