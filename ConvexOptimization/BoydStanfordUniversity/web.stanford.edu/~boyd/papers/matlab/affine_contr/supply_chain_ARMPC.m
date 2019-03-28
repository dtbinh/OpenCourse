% Affine-recourse MPC applied to supply chain example. 
% Comparison with affine controller and certainty-equivalent MPC. 

clear all;
randn('state',0);
rand('state',0);

% data generation
T = 10; 
x0 = 0;

s = 8;
h = 1;
b = 4;
g = 4;

mu = 3*rand(T,1)-1.5; 
tmp = randn(T)/4;
Sigma = tmp*tmp';

M = 500;                % number of training samples 
m = 500;                % number of testing samples

% Case I: Affine controller 
% -------------------------
[F, u0, x_tilde, u_tilde, Pxw, Puw] = supply_chain_AR(x0, s, b, h, g, mu, Sigma, M);
W = -exp(mu*ones(1,m) + sqrtm(Sigma)*randn(T,m));
xaf = Pxw*W + x_tilde*ones(1,m);
uaf = Puw*W + u_tilde*ones(1,m);
cost_affine = (sum(max(h*xaf(1:end-1,:),-b*xaf(1:end-1,:))) + max(g*xaf(end,:),-b*xaf(end,:)) + sum(abs(s*uaf)));


% Case II: Certainty-equivalent MPC
% ----------------------------------
cvx_quiet(true);
cost_mpc_ce = zeros(1,m); 
for iter = 1:m
    % at each t, replace future d-trajectory with mean of cond. distr.
    xce = zeros(T+1,1); 
    xce(1) = x0;
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
            x_0 = xce(t)*ones(T-t+2,1); 
            x = -G*d_bar + H*u + x_0; 
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


% Case III: Affine-recourse MPC
% ----------------------------
cvx_quiet(true);
cost_mpc_ar = zeros(1,m); 
for iter = 1:m    
    % at each t, sample from conditional distribution of d|d(0)...d(t-1)
    xar = zeros(T+1,1); 
    xar(1) = x0;
    uar = zeros(T,1);
    for t=1:T 
        disp(['AR-MPC *** iter = ' num2str(iter) ' t = ' num2str(t)]); 
        if t>1 
            v = mu(t:T) + Sigma(t:T,1:t-1)*(Sigma(1:t-1,1:t-1)\(log(-W(1:t-1,iter)) - mu(1:t-1)));
            V = Sigma(t:T,t:T) - Sigma(t:T,1:t-1)*(Sigma(1:t-1,1:t-1)\Sigma(1:t-1,t:T));
        else 
            v = mu;
            V = Sigma;
        end         
        [F_ar, u0_ar, x_tilde, u_tilde, Pxw, Puw, cvx_status] = supply_chain_AR(xar(t), s, b, h, g, v, V, M);
        if isempty(strfind(cvx_status, 'Solved')) 
            break
        end
        uar(t) = u0_ar(1) + F_ar(1,1)*xar(t); 
        xar(t+1) = xar(t) + uar(t) + W(t,iter);
    end  
    if isempty(strfind(cvx_status, 'Solved')) 
        cost_mpc_ar(iter) = NaN;
    else 
        cost_mpc_ar(iter) = sum(max(h*xar(1:end-1),-b*xar(1:end-1))) + max(g*xar(end),-b*xar(end)) + sum(abs(s*uar));
    end
end 

% Case IV: Compute prescient lower bound
% ---------------------------------------
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

disp(['affine: mean ' num2str(mean(cost_affine)) ', std '  num2str(std(cost_affine))]);
disp(['CE-MPC: mean ' num2str(mean(cost_mpc_ce(find(~isnan(cost_mpc_ce))))) ', std ' num2str(std(cost_mpc_ce(find(~isnan(cost_mpc_ce)))))]);
disp(['AR-MPC: mean ' num2str(mean(cost_mpc_ar(find(~isnan(cost_mpc_ar))))) ', std ' num2str(std(cost_mpc_ar(find(~isnan(cost_mpc_ar)))))]);
disp(['lwrbnd: mean ' num2str(mean(cost_lb)) ', std '  num2str(std(cost_lb))]);

figure;
specs = [50 450 0 50];
bins = linspace(50,800,100);
subplot(411)
set(gca, 'FontSize',18);
hist(cost_affine, bins);
hold on; 
line(mean(cost_affine)*ones(1,2), [0 100]);
axis(specs);
xlabel('affine');
subplot(412)
set(gca, 'FontSize',18);
hist(cost_mpc_ce, bins);
hold on; 
line(mean(cost_mpc_ce(find(~isnan(cost_mpc_ce))))*ones(1,2), [0 100]);
axis(specs);
xlabel('mpc-ce');
subplot(413)
set(gca, 'FontSize',18);
hist(cost_mpc_ar, bins);
hold on; 
line(mean(cost_mpc_ar(find(~isnan(cost_mpc_ar))))*ones(1,2), [0 100]);
axis(specs);
xlabel('mpc-ar');
subplot(414) 
set(gca, 'FontSize',18);
hist(cost_lb, bins);
hold on;
line(mean(cost_lb)*ones(1,2), [0 250]);
axis(specs);
xlabel('lower bound');
