close all
%load constr_tracking_sat_3

% Results
% -------
disp(['affine: mean ' num2str(mean(cost_affine)) ', std '  num2str(std(cost_affine))]);
disp(['affine (sat): mean ' num2str(mean(cost_sat_affine)) ', std '  num2str(std(cost_sat_affine))]);
disp(['nonlin: mean ' num2str(mean(cost_nonlin)) ', std '  num2str(std(cost_nonlin))]);
disp(['nonlin (sat): mean ' num2str(mean(cost_sat_nonlin)) ', std '  num2str(std(cost_sat_nonlin))]);
disp(['ce-mpc: mean ' num2str(mean(cost_mpc)) ', std '  num2str(std(cost_mpc))]);
disp(['lwrbnd: mean ' num2str(mean(cost_lb)) ', std '  num2str(std(cost_lb))]);

% Figures
% -------
figure 
idx = 7468; %round(rand*Mtest);
subplot(211)
set(gca, 'FontSize',14);
stairs(1:T, u_affine(:,idx),'--', 'LineWidth',2); 
hold on;
stairs(1:T, u_mpc(:,idx),'-.', 'LineWidth',2); 
stairs(1:T, u_nonlin(:,idx), 'LineWidth',2); 
stairs(1:T, u_lb(:,idx),'r', 'LineWidth',2); 
plot([1 T], u_max*ones(1,2),'k:')
plot([1 T], -u_max*ones(1,2),'k:')
ylabel('u') 
axis([1 T -1 1])

subplot(212)
set(gca, 'FontSize',14);
stairs(1:T, y_affine(:,idx),'--', 'LineWidth',2); 
hold on;
stairs(1:T, y_mpc(:,idx),'-.', 'LineWidth',2); 
stairs(1:T, y_nonlin(:,idx), 'LineWidth',2); 
stairs(1:T, y_lb(:,idx),'r', 'LineWidth',2); 
ylabel('y');
axis([1 T -2 2])
legend('affine', 'ce-mpc', 'nonlin', 'lwrbnd', 'Location', 'BestOutside', 'Orientation', 'horizontal'); 

print -depsc constr_tracking_sample_traj
return 

figure
specs = [0 30 0 650];
bins = linspace(0,35,100);
subplot(411)
% set(gca, 'FontSize',18);
hist(cost_affine, bins);
xlabel('affine');
hold on;
line(mean(cost_affine)*ones(1,2), [0 2000]);
axis(specs);

subplot(412)
% set(gca, 'FontSize',18);
hist(cost_mpc, bins);
xlabel('ce-mpc');
hold on;
line(mean(cost_mpc)*ones(1,2), [0 2000]);
axis(specs);

subplot(413)
% set(gca, 'FontSize',18);
hist(cost_nonlin, bins);
xlabel('nonlinear');
hold on;
line(mean(cost_nonlin)*ones(1,2), [0 2000]);
axis(specs);

subplot(414)
% set(gca, 'FontSize',18);
hist(cost_lb, bins);
xlabel('lower bound');
hold on;
line(mean(cost_lb)*ones(1,2), [0 2000]);
axis(specs);

print -depsc constr_tracking_err_hist

% figure;
% specs = [0 .8 0 15000];
% subplot(211)
% set(gca, 'FontSize',18);
% hist(abs(u_affine(:)), 75); 
% % line(mean(abs(u_affine(:)))*ones(1,2), [0 20000]);
% line(u_max*ones(1,2), [0 20000]);
% axis(specs);
% 
% subplot(212)
% set(gca, 'FontSize',18);
% hist(abs(u_nonlin(:)), 75); 
% % line(mean(abs(u_nonlin(:)))*ones(1,2), [0 20000]);
% line(u_max*ones(1,2), [0 20000]);
% axis(specs);
% 
% print -depsc constr_tracking_amp_hist

u_affine_short = u_affine(1:end-1,:); 
u_nonlin_short = u_nonlin(1:end-1,:); 

figure;
specs = [0 .8 0 15000];
subplot(211)
set(gca, 'FontSize',14);
hist(abs(u_affine_short(:)), 50); 
% line(mean(abs(u_affine(:)))*ones(1,2), [0 20000]);
line(u_max*ones(1,2), [0 20000]);
axis(specs);
xlabel('affine')

subplot(212)
set(gca, 'FontSize',14);
hist(abs(u_nonlin_short(:)), 50); 
% line(mean(abs(u_nonlin(:)))*ones(1,2), [0 20000]);
line(u_max*ones(1,2), [0 20000]);
axis(specs);
xlabel('nonlinear')

print -depsc constr_tracking_amp_hist_2
