% % Nonlinear Q-design for convex stochastic control 
% Section 6: Example (Tracking problem) 
% Consider a tracking problem with a scalar input, disturbance, sensor, and
% output signals. The tracking error, which is also the sensor signal, is
% given by:     y = Pyu*u - w
% Our objective is the mean-square tracking error E sum_{t=1}^T z_t^2. 
% The actuator input signal u_t must satisfy the (almost sure) constraint 
%               |u_t| <= u_max,     t=1,..., T. 
% The basis for \tilde Q will use levels \beta_1 = .5 and \beta_2 = 1. See
% section 5.3 and 6.1 for details. 

clear all;
randn('state',0);
rand('state',0);
cvx_quiet(true);

% data generation
% ---------------
T = 10;                         % time horizon
Pyw = -eye(T);
Pyu = zeros(T);
tau = .25;
for i=1:T 
    for j=1:i-1 
        Pyu(i,j) = exp(-tau*(i-j))*randn; 
    end
end
u_max = 0.5;                    % maximum actuator effort

% generate exogenous input sample trajectories
% --------------------------------------------
% variance of exogenous input vector
Sigma_w = zeros(T);             % variance of w
for i=1:T 
    for j=1:T
        Sigma_w(i,j) = exp(-((i-j)/3)^2); 
    end
end
L = chol(Sigma_w)'; 

Mtrain = 2000;                  % size of training set
Mtest = 10000;                  % size of validation set
wtrain = L*randn(T,Mtrain);     % validation set 
wtest = L*randn(T,Mtest);       % training set 

% standardized signals 
% --------------------
etrain = L\(Pyw*wtrain);                  
etest = L\(Pyw*wtest); 

% applying nonlinear bases
% ------------------------
betas = [0.5 1]; 
Qe_nonlin_train = cell(length(betas));
Qe_nonlin_test = cell(length(betas));
for i=1:length(betas)
    Qe_nonlin_train{i} =  max(min(etrain/betas(i),1),-1); 
end
for i=1:length(betas)
    Qe_nonlin_test{i} =  max(min(etest/betas(i),1),-1); 
end

% find nonlinear controller
% -------------------------
display('computing nonlinear controller ...'); 
cvx_begin
    variables F(T,T,length(betas))

    u = zeros(T,Mtrain);
    for i=1:length(betas) 
        u = u + F(:,:,i)*Qe_nonlin_train{i}; 
    end
    y = Pyw*wtrain + Pyu*u;

    % objective is square norm of tracking error 
    minimize (sum(square_pos(norms(y)))/Mtrain)
    
    % actuator effort constraint 
    abs(u) <= u_max; 

    % causality constraints
    for i=1:length(betas) 
        triu(F(:,:,i),1) == 0;
    end
cvx_end

% find optimal linear controller 
% ------------------------------
display('computing optimal linear controller ...'); 
cvx_begin
    variables Q_linear(T,T)

    u = Q_linear*etrain;
    y = Pyw*wtrain + Pyu*u;

    % objective is square norm of tracking error
    minimize (sum(square_pos(norms(y)))/Mtrain)
       
    % actuator effort constrains
    abs(u) <= u_max;

    % causality constraints
    triu(Q_linear,1) == 0;
cvx_end


% RESULTS ON VALIDATION SET 
% =========================
display('*** Testing controllers on the validation set ***'); 
display('-------------------------------------------------'); 

% test nonlinear controller
% -------------------------
display('testing nonlinear controller...'); 
u_nonlin = zeros(T,Mtest);
for i=1:length(betas)
    u_nonlin = u_nonlin + F(:,:,i)*Qe_nonlin_test{i};
end
y_nonlin = Pyu*u_nonlin + Pyw*wtest;
cost_nonlin = square_pos(norms(y_nonlin));

u_sat_nonlin = max(min(u_nonlin, u_max),-u_max);
y_sat_nonlin = Pyu*u_sat_nonlin + Pyw*wtest;
cost_sat_nonlin = square_pos(norms(y_sat_nonlin));

% test linear controller 
% ----------------------
display('testing optimal linear controller...'); 
u_linear = Q_linear*etest; 
y_linear = Pyu*u_linear + Pyw*wtest; 
cost_linear = square_pos(norms(y_linear));

u_sat_linear = max(min(u_linear, u_max),-u_max); 
y_sat_linear = Pyu*u_sat_linear + Pyw*wtest; 
cost_sat_linear = square_pos(norms(y_sat_linear)); 

% test CE-MPC
% -----------
display('testing CE-MPC (this will take a while!) ...');
cvx_quiet(true)
for m=1:Mtest 
    for t=1:T
        if t>1 
            w_hat = Sigma_w(t:T,1:t-1)*(Sigma_w(1:t-1,1:t-1)\wtest(1:t-1,m));
        else 
            w_hat = zeros(T,1);
        end        
        cvx_begin
            variable u_sub(T-t+1)
            y_sub = Pyu(t:end, t:end)*u_sub + Pyw(t:end,t:end)*w_hat; 
            minimize square_pos(norm(y_sub))
            abs(u_sub) <= u_max; 
        cvx_end
        u_mpc(t,m) = u_sub(1);
    end
end
y_mpc = Pyu*u_mpc + Pyw*wtest;
cost_mpc = square_pos(norms(y_mpc));

% compute prescient controller and bound
% --------------------------------------
display('testing prescient controller ...');
u_lb = zeros(T,Mtest);
y_lb = zeros(T,Mtest);
for m = 1:Mtest
    cvx_begin
        variable ulb(T)
        ylb = Pyu*ulb + Pyw*wtest(:,m);
        minimize (square_pos(norm(ylb)))
        abs(ulb) <= u_max
    cvx_end
    u_lb(:,m) = ulb;
    y_lb(:,m) = ylb;
end
cost_lb = square_pos(norms(y_lb));

% Results
% -------
disp('RESULTS ON THE VALIDATION SET: ');
disp('------------------------------'); 
disp(['linear: mean ' num2str(mean(cost_linear)) ', std '  num2str(std(cost_linear))]);
disp(['linear (sat): mean ' num2str(mean(cost_sat_linear)) ', std '  num2str(std(cost_sat_linear))]);
disp(['nonlin: mean ' num2str(mean(cost_nonlin)) ', std '  num2str(std(cost_nonlin))]);
disp(['nonlin (sat): mean ' num2str(mean(cost_sat_nonlin)) ', std '  num2str(std(cost_sat_nonlin))]);
disp(['ce-mpc: mean ' num2str(mean(cost_mpc)) ', std '  num2str(std(cost_mpc))]);
disp(['prescient: mean ' num2str(mean(cost_lb)) ', std '  num2str(std(cost_lb))]);

