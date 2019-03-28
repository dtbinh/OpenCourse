% This script generates all the figures
% for the paper:
% "Network Utility Maximization with Delivery Contracts"
% by N. Trichakis, A. Zymnis, and S. Boyd
% requires CVX to run

clear all; close all; clc;
cvx_quiet(true);
rand('state', 4);

%% --------- Problem Data -------------
% timesteps, links, flows
T = 10;
m = 3;
n = 3;

% routing matrix (is not time varying)
R = [1 0 1; 
     1 1 1; 
     0 1 1];

% capacity matrix
c_bar = [5;7;5]*ones(1, T); % mean capacities
delc= [1;3;1]*ones(1,T); % deviation +/- around mean
c = c_bar + delc.*(2*rand(m,T)-1);

% flow upper bounds
F_max = 4.5*ones(n,T);

% contracts
C1 = [ones(1,3) zeros(1,7); 
      zeros(1,5) ones(1,3) zeros(1,2)];
q1 = [12; 10];

C2 = [zeros(1,2) ones(1,4) zeros(1,4)];
q2 = 12;

C3 = [zeros(1,2) ones(1,8)];
q3 = 12;

%% ---------- The CVX Solver Call -----------
fprintf(1,'Solving with CVX...\n')
cvx_begin
    
    % the rate matrix F
    variable F(n, T)
    
    %the dual variables
    dual variable lambda
    dual variables mu1 mu2 mu3
    
    % the DNUM
    maximize    geomean(F(:))
    subject to  
                % link capacity constraints
                lambda: R*F <= c
                
                % contract constraints
                mu1: C1*F(1,:)' >= q1
                mu2: C2*F(2,:)' >= q2
                mu3: C3*F(3,:)' >= q3
                
                % positivity constraints
                F >= 0 
                F <= F_max
cvx_end
fprintf(1,'Done!\n\n')
Uopt = sum(log(F(:)));

% convert dual variables for sum(log(F))
lambda = lambda*(n*T)/cvx_optval;
mu1 = mu1*(n*T)/cvx_optval;
mu2 = mu2*(n*T)/cvx_optval;
mu3 = mu3*(n*T)/cvx_optval;

% create P = link_prices - subsidies
link_prices = R'*lambda;
subsidies = [C1'*mu1 C2'*mu2 C3'*mu3];
P = link_prices - subsidies';

%% ---------- Dual Decomposition Solution -----------
alpha = 0.01; %step size
Fdd = zeros(n,T);
lambda_dd = zeros(m,T);
mu1_dd = zeros(2,1); mu2_dd = 0; mu3_dd = 0;

fprintf(1,'Solving via Dual Decomposition...\n')
DDiters = 500;
for iter = 1:DDiters
    % create P = link_prices - subsidies
    link_prices = R'*lambda_dd;
    subsidies = [C1'*mu1_dd C2'*mu2_dd C3'*mu3_dd];
    P_dd = link_prices - subsidies';
    
    % update flows
    Fdd = min(1./P_dd,F_max);
    Fdd(find(P_dd<=0)) = F_max(find(P_dd<=0));
    L(iter) = sum(log(Fdd(:)));
    
    % update lambda
    cap_slack = c - R*Fdd;
    cap_viols(iter) = max(max(0,-cap_slack(:)));
    L(iter) = L(iter) + sum(sum(lambda_dd.*cap_slack));
    lambda_dd = max(lambda_dd-alpha*cap_slack,0);
    
    % update mu
    contr_slack1 = C1*Fdd(1,:)'-q1;
    contr_slack2 = C2*Fdd(2,:)'-q2;
    contr_slack3 = C3*Fdd(3,:)'-q3;
    contr_slacks = [contr_slack1; contr_slack2; contr_slack3];
    contr_viols(iter) = max(max(0,-contr_slacks));
    L(iter) = L(iter)+[mu1_dd; mu2_dd; mu3_dd]'*contr_slacks;
    mu1_dd = max(0,mu1_dd-alpha*contr_slack1);
    mu2_dd = max(0,mu2_dd-alpha*contr_slack2);
    mu3_dd = max(0,mu3_dd-alpha*contr_slack3);
end
fprintf(1,'Done!\n\n')

% compare the two solutions
diff_F = norm(F-Fdd,'fro');
diff_lambda = norm(lambda-lambda_dd,'fro')
diff_mu = norm([mu1;mu2;mu3]-[mu1_dd;mu2_dd;mu3_dd]);
fprintf(1,'||F - Fdd|| = %3.3f, ||lambda - lambda_dd|| = %3.3f, ||mu-mu_dd|| = %3.3f\n\n',...
    diff_F,diff_lambda,diff_mu);

%% ---------- MPC heuristic using CVX -----------
fprintf(1,'MPC heuristic using CVX...\n')
F_bar = [];
for t = 1:T
    fprintf(1,'MPC timestep %d of %d...\n',t,T);
    
    c_hat = [c(:,1:t) c_bar(:,t+1:T)];
    cvx_begin

        % the rate matrix F
        variable Fmpc(n, T)

        % the DNUM
        F2 = Fmpc(:,t:end);
        maximize    geomean(F2(:))
        subject to
            % fix rates prior to time t
            if t>1
                Fmpc(:,1:t-1) == F_bar(:,1:t-1);
            end
            
            % link capacity constraints
            R*Fmpc <= c_hat

            % contract constraints
            C1*Fmpc(1,:)' >= q1
            C2*Fmpc(2,:)' >= q2
            C3*Fmpc(3,:)' >= q3

            % positivity constraints
            Fmpc >= 0
            Fmpc <= F_max
    cvx_end
    
    % update rates
    F_bar = Fmpc;
end        
fprintf(1,'Done!\n\n')

%% --------------- Generate Figures -----------

% convergence plot
figure
plot(1:DDiters, Uopt*ones(1,DDiters),'k--','linewidth',2); hold on;
set(gca,'FontSize',12);
plot(L,'linewidth',2);
xlabel('iter');
ylabel('upperbound');
axis([0 DDiters 0 ceil(max(L))])
box off;
print -depsc ubound.eps

% violation plots
figure
cap_viols(find(cap_viols==0))=1e-4;
subplot(2,1,1);
semilogy(1:DDiters, cap_viols, 'linewidth',2);
set(gca,'FontSize',12);
ylabel('capacity');
axis([1 DDiters 1e-2 1e2]);
box off

contr_viols(find(contr_viols==0))=1e-4;
subplot(2,1,2);
semilogy(1:DDiters, contr_viols, 'linewidth',2);
set(gca,'FontSize',12);
xlabel('iter');
ylabel('contract');
axis([1 DDiters 1e-2 1e2]);
box off
print -depsc violations.eps

% plot optimal rates
figure; colors = get(gca,'ColorOrder');
v_max = ceil(max(F(:)));
subplot(3,1,1); hold on
fill([0 3.5 3.5 0],[0 0 v_max v_max],0.95*[1 1 1],'EdgeColor','none');
fill([5.5 8.5 8.5 5.5],[0 0 v_max v_max],0.95*[1 1 1],'EdgeColor','none');
stairs(0.5:1:T+.5, [F(1,:) F(1,end)],'-','Color',colors(1,:),'LineWidth',2);
set(gca,'FontSize',12);
ylabel('f1');
axis([.5 T+.5 0 v_max]);
box off;

subplot(3,1,2); hold on
fill([2.5 6.5 6.5 2.5],[0 0 v_max v_max],0.95*[1 1 1],'EdgeColor','none');
stairs(0.5:1:T+.5, [F(2,:) F(2,end)],'-','Color',colors(2,:),'LineWidth',2);
set(gca,'FontSize',12);
ylabel('f2');
axis([.5 T+.5 0 v_max]);
box off;

subplot(3,1,3); hold on
fill([2.5 10.5 10.5 2.5],[0 0 v_max v_max],0.95*[1 1 1],'EdgeColor','none');
stairs(0.5:1:T+.5, [F(3,:) F(3,end)],'-','Color',colors(3,:),'LineWidth',2);
set(gca,'FontSize',12);
ylabel('f3');
axis([.5 T+.5 0 v_max]);
xlabel('time step');
box off;
print -depsc rate.eps

% plot link traffic
traffic = R*F; v_max = ceil(max(c(:)));
figure; subplot(3,1,1); hold on;
stairs([0.5:1:T+.5],[c(1,:) c(1,end)],'k-.','LineWidth',2);
stairs([0.5:1:T+.5],[traffic(1,:) traffic(1,end)],'LineWidth',2);
ylabel('link1');
axis([0.5 T+0.5 0 v_max]);
set(gca,'FontSize',12);
box off;

subplot(3,1,2); hold on;
stairs([0.5:1:T+.5],[c(2,:) c(2,end)],'k-.','LineWidth',2);
stairs([0.5:1:T+.5],[traffic(2,:) traffic(2,end)],'LineWidth',2);
ylabel('link2');
axis([0.5 T+0.5 0 v_max]);
set(gca,'FontSize',12);
box off;

subplot(3,1,3); hold on;
stairs([0.5:1:T+.5],[c(3,:) c(3,end)],'k-.','LineWidth',2);
stairs([0.5:1:T+.5],[traffic(3,:) traffic(3,end)],'LineWidth',2);
ylabel('link3');
axis([0.5 T+0.5 0 v_max]);
set(gca,'FontSize',12);
box off;
xlabel('time step');
print -depsc linktraffic.eps

% plot cumulative rates
Fcum = cumsum(F,2); colors = get(gca,'ColorOrder'); 
v_max = ceil(max(Fcum(:)));
figure; subplot(3,1,1); hold on
plot([0 3.5],[0 q1(1)],'o--','Color',colors(1,:),'LineWidth',2);
plot([5.5 8.5],[Fcum(1,5) Fcum(1,5)+q1(2)],'o--','Color',colors(1,:),'LineWidth',2);
stairs([0.5:1:T+.5], [Fcum(1,:) Fcum(1,end)],'-','Color',colors(1,:),'LineWidth',2);
set(gca,'FontSize',12);
ylabel('cumrate1');
axis([.5 T+.5 0 v_max]);
box off;

subplot(3,1,2); hold on
plot([2.5 6.5],[Fcum(2,2) Fcum(2,2)+q2],'o--','Color',colors(2,:),'LineWidth',2);
stairs([0.5:1:T+.5], [Fcum(2,:) Fcum(2,end)],'-','Color',colors(2,:),'LineWidth',2);
set(gca,'FontSize',12);
ylabel('cumrate2');
axis([.5 T+.5 0 v_max]);
box off;

subplot(3,1,3); hold on
plot([2.5 10.5],[Fcum(3,2) Fcum(3,2)+q3],'o--','Color',colors(3,:),'LineWidth',2);
stairs([0.5:1:T+.5], [Fcum(3,:) Fcum(3,end)],'-','Color',colors(3,:),'LineWidth',2);
set(gca,'FontSize',12);
ylabel('cumrate3');
axis([.5 T+.5 0 v_max]);
box off;
xlabel('time step');
print -depsc cumrate.eps

% plot flow prices
figure; colors = get(gca,'ColorOrder'); 
v_max = ceil(max(10*P(:)))/10;
v_min = min(0,floor(min(10*P(:)))/10);
subplot(3,1,1); hold on
fill([0 3.5 3.5 0],[v_min v_min v_max v_max],0.95*[1 1 1],'EdgeColor','none');
fill([5.5 8.5 8.5 5.5],[v_min v_min v_max v_max],0.95*[1 1 1],'EdgeColor','none');
stairs(0.5:1:T+.5, [P(1,:) P(1,end)],'-','Color',colors(1,:),'LineWidth',2);
set(gca,'FontSize',12);
ylabel('price1');
axis([.5 T+.5 v_min v_max]);
box off;

subplot(3,1,2); hold on
fill([2.5 6.5 6.5 2.5],[v_min v_min v_max v_max],0.95*[1 1 1],'EdgeColor','none');
stairs(0.5:1:T+.5, [P(2,:) P(2,end)],'-','Color',colors(2,:),'LineWidth',2);
set(gca,'FontSize',12);
ylabel('price2');
axis([.5 T+.5 v_min v_max]);
box off;

subplot(3,1,3); hold on
fill([2.5 10.5 10.5 2.5],[v_min v_min v_max v_max],0.95*[1 1 1],'EdgeColor','none');
stairs(0.5:1:T+.5, [P(3,:) P(3,end)],'-','Color',colors(3,:),'LineWidth',2);
set(gca,'FontSize',12);
ylabel('price3');
axis([.5 T+.5 v_min v_max]);
xlabel('time step');
box off;
print -depsc prices.eps

% plot link prices
v_max = ceil(10*max(lambda(:)))/10;
figure; subplot(3,1,1); hold on;
stairs([0.5:1:T+.5],[lambda(1,:) lambda(1,end)],'b-','LineWidth',2);
ylabel('link1');
axis([0.5 T+0.5 0 v_max]);
set(gca,'FontSize',12);
box off;

subplot(3,1,2); hold on;
stairs([0.5:1:T+.5],[lambda(2,:) lambda(2,end)],'b-','LineWidth',2);
ylabel('link2');
axis([0.5 T+0.5 0 v_max]);
set(gca,'FontSize',12);
box off;

subplot(3,1,3); hold on;
stairs([0.5:1:T+.5],[lambda(3,:) lambda(3,end)],'b-','LineWidth',2);
ylabel('link3');
axis([0.5 T+0.5 0 v_max]);
set(gca,'FontSize',12);
xlabel('time step');
box off;
print -depsc lambda.eps

% plot MPC and prescient comparison
figure; colors = get(gca,'ColorOrder');
v_max = ceil(max(Fmpc(:)));
subplot(3,1,1); hold on
fill([0 3.5 3.5 0],[0 0 v_max v_max],0.95*[1 1 1],'EdgeColor','none');
fill([5.5 8.5 8.5 5.5],[0 0 v_max v_max],0.95*[1 1 1],'EdgeColor','none');
stairs(0.5:1:T+.5, [F(1,:) F(1,end)],'--','Color',colors(1,:),'LineWidth',2);
stairs(0.5:1:T+.5, [Fmpc(1,:) Fmpc(1,end)],'-','Color',colors(1,:),'LineWidth',2);
set(gca,'FontSize',12);
ylabel('f1');
axis([.5 T+.5 0 v_max]);
box off;

subplot(3,1,2); hold on
fill([2.5 6.5 6.5 2.5],[0 0 v_max v_max],0.95*[1 1 1],'EdgeColor','none');
stairs(0.5:1:T+.5, [F(2,:) F(2,end)],'--','Color',colors(2,:),'LineWidth',2);
stairs(0.5:1:T+.5, [Fmpc(2,:) Fmpc(2,end)],'-','Color',colors(2,:),'LineWidth',2);
set(gca,'FontSize',12);
ylabel('f2');
axis([.5 T+.5 0 v_max]);
box off;

subplot(3,1,3); hold on
fill([2.5 10.5 10.5 2.5],[0 0 v_max v_max],0.95*[1 1 1],'EdgeColor','none');
stairs(0.5:1:T+.5, [F(3,:) F(3,end)],'--','Color',colors(3,:),'LineWidth',2);
stairs(0.5:1:T+.5, [Fmpc(3,:) Fmpc(3,end)],'-','Color',colors(3,:),'LineWidth',2);
set(gca,'FontSize',12);
ylabel('f3');
axis([.5 T+.5 0 v_max]);
xlabel('time step');
box off;
print -depsc rate_mpc.eps
    
    