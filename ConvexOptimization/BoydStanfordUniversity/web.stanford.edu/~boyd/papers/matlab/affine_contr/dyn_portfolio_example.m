% Dynamic portfolio optimization: income investing 
% Example from "Design of Affine Controllers via Convex Optimization",
% Joelle Skaf and Stephen Boyd.
% www.stanford.edu/~boyd/papers/affine_contr.html
% Requires CVX package to run.

rand('state',0);
rand('state',0);

% data generation 
% ---------------
T = 10;                 % number of periods
n = 4;                  % number of assets 
kappa = 0.1;            % transaction cost rate
rho = .1; 
r_min = 0.10;

% return r(t) ~ exp(N(mu, Sigma))-1
beta = 0.3+0.7*rand(n,1);           % beta's are unif on [0.3,1]
beta(1) = 0;                        % beta for risk-free asset
sigma_m = 0.1;                      % market std dev
sigmas = 0.2*rand(n,1);             % firm-specific std dev unif on [0, 20]
sigmas(1) = 0;                      % std dev for risk-free asset 
Sigma = sigma_m^2*beta*beta' + diag(sigmas.^2);

mu_rf = 0.02;                       % risk-free return
SR = 0.4;                           % sharpe ratio
mu = mu_rf + SR*sqrt(diag(Sigma));  % expected return is 0.4*std dev

[mu, idx] = sort(mu, 'ascend');
Sigma = sigma_m^2*beta(idx)*beta(idx)' + diag(sigmas(idx).^2);

rbar = exp(mu+0.5*diag(Sigma)); 
Sigma_bar = (rbar*rbar').*(exp(Sigma)-1);
rbar = rbar - 1; 

% target portfolio 
% ----------------
xtar = [1;0;1;1]/3;

% affine controller 
% -----------------
I = eye(n);
A = I + diag(rbar);
B = I + diag(rbar); 
x_0 = xtar; 

x0 = zeros(n*(T+1),1);
tmp = x_0;
for i=0:T
    x0(n*i+1:n*(i+1)) = tmp;
    tmp = A*tmp;
end

G = zeros(n*(T+1),n*T); 
tmp = eye(n); 
G(n+1:2*n,1:n) = tmp; 
for i=2:T
    tmp = A*tmp; 
    G(n*i+1:n*(i+1), 1:n*i) = [tmp G(n*(i-1)+1:n*i, 1:n*(i-1))]; 
end
   
H = zeros(n*(T+1),n*T); 
tmp = B; 
H(n+1:2*n,1:n) = tmp; 
for i=2:T
    tmp = A*tmp; 
    H(n*i+1:n*(i+1), 1:n*i) = [tmp H(n*(i-1)+1:n*i, 1:n*(i-1))]; 
end

S = zeros(n*T); 
for i=0:T-1
    S(n*i+1:n*i+n, n*i+1:n*i+n) = diag(xtar)*sqrtm(Sigma_bar); 
end

% generate M random samples of r
%M = 1000;  % used in paper.
M = 100;
W = [];
for t=0:T-1
    R = exp(mu*ones(1,M) + sqrtm(Sigma)*randn(n,M))-1;         
    W = [W; (R - rbar*ones(1,M)).*(xtar*ones(1,M))];
end

I = eye(n*(T+1));
cvx_begin
    variable Q(n*T, n*(T+1))
    variable r(n*T) 
    
    Pxw = (I+H*Q)*G;
    Puw = Q*G;
    x_tilde = (I + H*Q)*x0 + H*r;
    u_tilde = Q*x0 + r;
    
    x = Pxw*W + x_tilde*ones(1,M); 
    u = Puw*W + u_tilde*ones(1,M);   

    obj = 0; 
    for t=0:T-1
        ut = u(n*t+1:n*(t+1), :);
        obj = obj -square_pos (r_min+sum(ut) + kappa*norms(ut,1));
    end
    maximize (sum(obj)/(M*T))
    
    % constraints 
    Z = Pxw*S;       
    for t=0:T
        Zt = Z(n*t+1:n*(t+1), :);
        square_pos(norm(Zt,'fro')) + square_pos(norm(x_tilde(n*t+1:n*(t+1), :) -xtar)) <= n*rho^2
        
        sum(x_tilde(n*t+1:n*(t+1), :)) == 1;
    end
    
    % Q (n,n) block lower triangular
    for i=0:T-1 
        Q(n*i+1:n*(i+1), n*(i+1)+1:end) == 0
    end    
cvx_end

I = eye(n*T);
F = (I+Q*H)\Q;
u0 = (I+Q*H)\r;

% test affine controller 
% ----------------------
% generate m random samples of r
m = 5000;
Wt = [];
Rs = []; 
for t=0:T-1
    R = exp(mu*ones(1,m) + sqrtm(Sigma)*randn(n,m))-1;         
    Rs = [Rs; R]; 
    Wt = [Wt; (R - rbar*ones(1,m)).*(xtar*ones(1,m))];    
end

% propagating using linearized dynamics 
x_linear = Pxw*Wt + x_tilde*ones(1,m);
u_linear  = Puw*Wt + u_tilde*ones(1,m);
util_linear = 0;
cash_linear = zeros(T,m);
for t=0:T-1
    ut = u_linear(n*t+1:n*(t+1), :);
    util_linear = util_linear -square_pos (r_min+sum(ut) + kappa*norms(ut,1));
    cash_linear(t+1,:) = -sum(ut) - kappa*norms(ut,1);
end
% util_linear = util_linear/T;
% cash_linear = -sum(u_linear) - kappa*norms(u_linear,1);

% propagating using nonlinear dynamics
x_nonlin = zeros(n*(T+1),m); 
u_nonlin = zeros(n*T, m); 
cash_nonlin = zeros(T,m);
for iter =1:m
    r = Rs(:,iter); 
    x = x_0; 
    xall = zeros(n,T+1);
    xall(:,1) = x_0; 
    uall = zeros(n,T);     
    for i=1:T
        tmp = xall(:,1:i); 
        tmp = tmp(:); 
        uall(:,i) = F(n*(i-1)+1:n*i, 1:n*i)*tmp + u0(n*(i-1)+1:n*i);
        xall(:,i+1) = (1+r(n*(i-1)+1:n*i)).*(xall(:,i) + uall(:, i));
    end
    x_nonlin(:, iter) = xall(:);
    u_nonlin(:, iter) = uall(:);
end
util_nonlin = 0;
for t=0:T-1
    ut = u_nonlin(n*t+1:n*(t+1), :);
    util_nonlin = util_nonlin -square_pos (r_min+sum(ut) + kappa*norms(ut,1));
    cash_nonlin(t+1,:) = -sum(ut) - kappa*norms(ut,1);
end
% util_nonlin = util_nonlin/T;
% cash_nonlin = -sum(u_nonlin) - kappa*norms(u_nonlin,1);

% test greedy nonlinear controller 
% --------------------------------
x2 = zeros(n*(T+1),m); 
u2 = zeros(n*T, m); 
cash2 = zeros(T,m);
D = diag(diag(Sigma_bar)) + diag((1+rbar).^2);
for iter = 1:m 
%     disp(iter);
    r = reshape(Rs(:,iter),n,T); 
    x = x0;
    xall = zeros(n,T+1);
    xall(:,1) = x_0; 
    uall = zeros(n,T);    
    for i=1:T 
        xprev = xall(:,i); 
        z = xtar./(1+rbar) - xprev;
        uall(:,i) = z;
        xall(:,i+1) = (1+r(:,i)).*(xprev + z);
    end
    x2(:,iter) = xall(:);
    u2(:,iter) = uall(:);
end
util2 = 0;
for t=0:T-1
    ut = u2(n*t+1:n*(t+1), :);
    util2 = util2 -square_pos (r_min+sum(ut) + kappa*norms(ut,1));
    cash2(t+1,:) = -sum(ut) - kappa*norms(ut,1);
end
% util2 = util2/T;
% cash2 = -sum(u2) - kappa*norms(u2,1);

% plots
% -----
figure;
specs = [-0.15 0 0 700];
bins = linspace(-0.16,0,100);
subplot(311); 
set(gca, 'FontSize',18);
hist(util_linear,bins);
hold on; 
line(mean(util_linear)*ones(1,2), [0 1000]);
axis(specs);
title('linearized dynamics')
subplot(312);
set(gca, 'FontSize',18);
hist(util_nonlin,bins);
hold on;
line(mean(util_nonlin)*ones(1,2), [0 1000]);
axis(specs);
title('nonlinear dynamics')
subplot(313);
set(gca, 'FontSize',18);
hist(util2,bins);
hold on;
line(mean(util2)*ones(1,2), [0 1000]);
axis(specs);
title('greedy controller')

figure;
specs = [-.2 .4 0 2500];
bins = linspace(-0.3,.5, 500);
subplot(311); 
set(gca, 'FontSize',18);
hist(cash_linear(:),bins);
hold on;
axis(specs);
title('linear dynamics')
subplot(312); 
set(gca, 'FontSize',18);
hist(cash_nonlin(:),bins);
hold on;
axis(specs);
title('nonlinear dynamics')
subplot(313);
set(gca, 'FontSize',18);
hist(cash2(:),bins);
hold on;
axis(specs);
title('greedy controller')
