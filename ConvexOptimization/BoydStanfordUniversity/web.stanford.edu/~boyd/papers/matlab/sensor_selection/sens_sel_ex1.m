% Generates the numerical example in the paper
% Sensor Selection via Convex Optimization
% www.stanford.edu/~boyd/papers/sensor_selection.html
%
% May 2008 Siddharth Joshi & Stephen Boyd
clear all;
m = 100; % number of sensors
n = 20; % dimension of x to be estimated
ks = 20:1:40; % target number of sensors to use

randn('state', 0);
A = randn(m,n)/sqrt(n);

fprintf('\nSensor selection problem:\nNumber of sensors: m = %d\nNumber of parameters: n = %d\n', m, n);
L = []; Utilde = []; L_loc = [];  L_loc2 = []; threshold = 0.4;
i = 0;
for k = ks
    fprintf('\n\nSelecting k = %d sensors...\n', k);
    i = i + 1;
    [zhat L(i) zast Utilde(i)] = sens_sel_approxnt(A, k);
    [z_loc L_loc(i)] = sens_sel_loc(A, zhat);
    [z_loc2 L_loc2(i)] = sens_sel_locr(A, k, zast, threshold);
    fprintf('\nUtilde: %.3e, L: %.3e,  L_loc: %.3e, L_loc2: %.3e\n', Utilde(i), L(i), L_loc(i), L_loc2(i));
end
delta = Utilde - L; delta_loc = Utilde - L_loc; delta_loc2 =  Utilde - L_loc2;


figure;
subplot(2,1,1); hold on;
set(gca,'FontName','times', 'FontSize', 16);
xlabel('k'); ylabel('bounds'); 
plot(ks, L, 'b-', 'LineWidth', 1);
plot(ks, L_loc, 'k--', 'LineWidth', 1);
plot(ks, L_loc2, 'm-.', 'LineWidth', 1);
plot(ks, Utilde, 'r:', 'LineWidth', 1);
hold off

subplot(2,1,2); hold on;
set(gca,'FontName','times', 'FontSize', 16);
xlabel('k'); ylabel('gaps');
plot(ks, delta, 'b-', 'LineWidth', 1);
plot(ks, delta_loc, 'k--', 'LineWidth', 1);
plot(ks, delta_loc, 'm-.', 'LineWidth', 1);
hold off
print -deps ex1.eps

figure; hold on;
set(gca,'FontName','times', 'FontSize', 16);
xlabel('k'); ylabel('radiiratio');
plot(ks, exp(delta/(2*n)), 'b-', 'LineWidth', 1);
plot(ks, exp(delta_loc/(2*n)), 'k--', 'LineWidth', 1);
plot(ks, exp(delta_loc2/(2*n)), 'm-.', 'LineWidth', 1);
hold off
print -deps ex1r.eps
