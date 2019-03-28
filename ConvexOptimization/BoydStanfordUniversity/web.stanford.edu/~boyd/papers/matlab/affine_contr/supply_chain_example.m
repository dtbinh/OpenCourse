% Linear supply chain example: x(t+1) = x(t) + u(t) - d(t).
% x, u and d are scalar.
% Find affine controller that minimizes sum(Phi(x(t)).
% here d ~ LN(mu, Sigma) (here d is whole trajectory) 
%
% Example from "Design of Affine Controllers via Convex Optimization",
% Joelle Skaf and Stephen Boyd.
% www.stanford.edu/~boyd/papers/affine_contr.html
% Requires CVX package to run.

clear all;
randn('state',0);
rand('state',0);

% data generation
T = 10; 

G = ones(T+1);
G = tril(G) - diag(diag(G)); 
G = G(:,1:end-1); 
H = G; 
x_0 = 0;
x0 = x_0*ones(T+1,1); 

s = 8;
h = 1;
b = 4;
g = 4;

mu = 3*rand(T,1)-1.5; 
tmp = randn(T)/4;
Sigma = tmp*tmp';


% Case I: optimal affine controller
% ---------------------------------
% number of sample paths for w 
% M = 5000; % used in paper; takes many minutes to run!
M = 500;
W = -exp(mu*ones(1,M) + sqrtm(Sigma)*randn(T,M));
I = eye(T+1); 
cvx_begin
    variable Q(T,T+1) 
    variable r(T) 
       
    Pxw = (I+H*Q)*G;
    Puw = Q*G;
    x_tilde = (I + H*Q)*x0 + H*r;
    u_tilde = Q*x0 + r;
    
    x = Pxw*W + x_tilde*ones(1,M);
    u = Puw*W + u_tilde*ones(1,M);
    minimize(sum(sum(max(h*x(1:end-1,:),-b*x(1:end-1,:))) + max(g*x(end,:),-b*x(end,:)) + sum(abs(s*u))))
    Qupper = [Q;zeros(1,T+1)];
    Qupper = triu(Qupper) - diag(diag(Qupper));
    Qupper == 0
cvx_end

I = eye(T);
F = (I+Q*H)\Q;
u0 = (I+Q*H)\r;

% testing the control laws obtained on new instances
%M = 5000; % used in paper
M = 200;
W = -exp(mu*ones(1,M) + sqrtm(Sigma)*randn(T,M));
% Case I: affine controller
x = Pxw*W + x_tilde*ones(1,M);
u = Puw*W + u_tilde*ones(1,M);
cost_affine = (sum(max(h*x(1:end-1,:),-b*x(1:end-1,:))) + max(g*x(end,:),-b*x(end,:)) + sum(abs(s*u)));

% Case II: naive greedy controller
% --------------------------------
cost_naive = zeros(1,M);
d_mean = exp(mu+.5*diag(Sigma));
for iter =1:M
    % simple controller  
    xn = [0; d_mean + W(:,iter)]; 
    un = d_mean - xn(1:end-1);
    cost_naive(iter) = sum(max(h*xn(1:end-1),-b*xn(1:end-1))) + max(g*xn(end),-b*xn(end)) + sum(abs(s*un));
end

% Case III: greedy controller
% ---------------------------
cost_greedy = zeros(1, M);
for iter =1:M
    % simple controller 
    % d_bar(t) is the expected value of d(t) given d(0),...,d(t-1)
    d_bar = zeros(T, 1); 
    d_bar(1) = exp(mu(1) + 0.5*Sigma(1,1));
    for t=2:T 
        v = mu(t) + Sigma(t,1:t-1)*(Sigma(1:t-1,1:t-1)\(log(-W(1:t-1,iter)) - mu(1:t-1)));
        V = Sigma(t,t) - Sigma(t,1:t-1)*(Sigma(1:t-1,1:t-1)\Sigma(1:t-1,t));
        d_bar(t) = exp(v+0.5*V);
    end    
    xs = [0; d_bar + W(:,iter)]; 
    us = d_bar - xs(1:end-1);
    cost_greedy(iter) = sum(max(h*xs(1:end-1),-b*xs(1:end-1))) + max(g*xs(end),-b*xs(end)) + sum(abs(s*us));
end

% Case IV: CE-MPC
% ---------------
cvx_quiet(true);
cost_mpc_ce = zeros(1,M); 
for iter = 1:M
    % at each t, replace future d-trajectory with mean of cond. distr.
    xce = zeros(T+1,1); 
    xce(1) = x_0;
    uce = zeros(T,1);
    for t=1:T 
        disp(['CE-MPC *** iter = ' num2str(iter) ' t = ' num2str(t)]); 
        if t>1 
            v = mu(t:T) + Sigma(t:T,1:t-1)*(Sigma(1:t-1,1:t-1)\(log(-W(1:t-1,iter)) - mu(1:t-1)));
            V = Sigma(t:T,t:T) - Sigma(t:T,1:t-1)*(Sigma(1:t-1,1:t-1)\Sigma(1:t-1,t:T));
        else 
            v = mu;
            V = Sigma;
        end
        d_bar = exp(v + 0.5*diag(V));       % mean of conditional distr.
        cvx_begin 
            variable u(T-t+1) 
            G = ones(T-t+2);
            G = tril(G) - diag(diag(G)); 
            G = G(:,1:end-1); 
            H = G; 
            xinit = xce(t)*ones(T-t+2,1); 
            x = -G*d_bar + H*u + xinit; 
            minimize sum(max(h*x(1:end-1),-b*x(1:end-1))) + max(g*x(end),-b*x(end)) + sum(abs(s*u))
        cvx_end 
        if isempty(strfind(cvx_status, 'Solved')) 
            break
        end
        uce(t) = u(1);
        xce(t+1) = xce(t) + uce(t) + W(t,iter);
    end  
    if isempty(strfind(cvx_status, 'Solved')) 
        cost_mpc_ce(iter) = NaN;
    else 
        cost_mpc_ce(iter) = sum(max(h*xce(1:end-1),-b*xce(1:end-1))) + max(g*xce(end),-b*xce(end)) + sum(abs(s*uce));
    end
end 

% Case V: Compute prescient lower bound 
% -------------------------------------
G = ones(T+1);
G = tril(G) - diag(diag(G)); 
G = G(:,1:end-1); 
H = G; 
cost_lb = zeros(1,M); 
for iter = 1:M
    disp(['Lower Bound *** iter = ' num2str(iter) ' t = ' num2str(t)]); 
    cvx_begin
        variable ulb(T) 
        xlb = G*W(:,iter) + H*ulb;
        minimize sum(max(h*xlb(1:end-1),-b*xlb(1:end-1))) + max(g*xlb(end),-b*xlb(end)) + sum(abs(s*ulb))
    cvx_end
    if isempty(strfind(cvx_status, 'Solved')) 
        cost_lb(iter) = NaN;
    else 
        cost_lb(iter) = sum(max(h*xlb(1:end-1),-b*xlb(1:end-1))) + max(g*xlb(end),-b*xlb(end)) + sum(abs(s*ulb));
    end
end

% Results
% -------
disp(['naive: mean ' num2str(mean(cost_naive)) ', std ' num2str(std(cost_naive))]);
disp(['greedy: mean ' num2str(mean(cost_greedy)) ', std ' num2str(std(cost_greedy))]);
disp(['CE-MPC: mean ' num2str(mean(cost_mpc_ce)) ', std ' num2str(std(cost_mpc_ce))]);
disp(['affine: mean ' num2str(mean(cost_affine)) ', std '  num2str(std(cost_affine))]);
disp(['lwrbnd: mean ' num2str(mean(cost_lb)) ', std '  num2str(std(cost_lb))]);

% Plots
% -----
figure;
specs = [0 700 0 100];
bins = linspace(50,800,100);
subplot(511) 
set(gca, 'FontSize',18);
hist(cost_naive, bins);
xlabel('naive greedy');
hold on;
line(mean(cost_naive)*ones(1,2), [0 250]);
axis(specs);
subplot(512) 
set(gca, 'FontSize',18);
hist(cost_greedy, bins);
xlabel('greedy');
hold on;
line(mean(cost_greedy)*ones(1,2), [0 250]);
axis(specs);
subplot(513) 
set(gca, 'FontSize',18);
hist(cost_mpc_ce, bins);
xlabel('CE-MPC');
hold on;
line(mean(cost_mpc_ce)*ones(1,2), [0 250]);
axis(specs);
subplot(514) 
set(gca, 'FontSize',18);
hist(cost_affine, bins);
xlabel('affine');
hold on;
line(mean(cost_affine)*ones(1,2), [0 250]);
axis(specs);
subplot(515) 
set(gca, 'FontSize',18);
hist(cost_lb, bins);
xlabel('lower bound');
hold on;
line(mean(cost_lb)*ones(1,2), [0 250]);
axis(specs);
