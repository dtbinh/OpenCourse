close all;

% PARAMETERS (mecrow)
addpath('mecrow_model')
[R, M, C, Np, f, g, Vmax] = mecrow_model();
motor = sw_reluc(R, M, C, f, g, Vmax);
F_tilde = [0, 1000, 2000, 3500, 6000];

% PRINT PLOTS OF THE PWA APPROXIMATIONS
if 0
  motor.plot_pwa_3d();
end

% (RE)-GENERATE PWA MODELS
if 0
  theta = linspace(0, 2*pi/(3*Np), 40);
  motor = motor.find_abcd(F_tilde, theta);
  a = motor.a;
  b = motor.b;
  c = motor.c;
  d = motor.d;
  T = motor.T;
  N = motor.N;
else
  motor.a = a;
  motor.b = b;
  motor.c = c;
  motor.d = d;
  motor.T = T;
  motor.N = N;
  motor.F_tilde = F_tilde;
  motor.theta = theta;
end

% OPTIMIZE FOR A SPECIFIC OPERATING POINT
mu = 5;
%motor = motor.find_D();  % try out 'improved' formulation
if 1
  omega = 1000 * 2*pi/60;
  tau_des = 10;
  [i, v, tau, z, q, F, psi, phi, lambda, time] = motor.optimize(mu, tau_des, omega);

  [theta3, i3] = motor.onephase2three(theta, i);
  [theta3, tau3] = motor.onephase2three(theta, tau);
  [theta3, v3] = motor.onephase2three(theta, v);

  close all;
  figure;
  subplot(311)
  plot(theta3, i3)
  ylabel('i')
  subplot(312)
  plot(theta3, v3)
  ylabel('v')
  subplot(313)
  plot(theta3, tau3, theta3, sum(tau3))
  ylabel('tau')
  xlabel('th')
  print -depsc sim_opt.eps
end

% CHECK TIMING FOR RANDOM OPERATING POINTS
if 0
  N_iter = 3;
  rng(0)
  tau_ub = 15;
  tau_lb = -15;
  omega_ub = 4000 * 2*pi/60;
  omega_lb = 0;
  tau_vec = (tau_ub - tau_lb)*rand(N_iter,1) + tau_lb;
  omega_vec = (omega_ub - omega_lb)*rand(N_iter,1) + omega_lb;
  for iter = 1:N_iter
    [~, ~, ~, ~, ~, ~, ~, ~, ~, time(iter)] = ... 
         motor.optimize(mu, tau_vec(iter), omega_vec(iter));
  end
  fprintf('average time:  %f\n\n', mean(time));
end

