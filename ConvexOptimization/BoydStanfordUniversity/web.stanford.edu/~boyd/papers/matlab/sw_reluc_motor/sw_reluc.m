classdef sw_reluc

properties
  R, M, C, f, g, Vmax, l, m, n, Q, D, a, b, c, d, F_tilde, T, N, theta
end

methods
  function s = sw_reluc(R, M, C, f, g, Vmax)
    s.R = R;
    s.M = M;
    s.C = C;
    s.f = f;
    s.g = g;
    s.Vmax = Vmax;
    [s.l, s.n] = size(C);
    s.m = size(M, 2);
  end

  function s = find_abcd(s, F_tilde, theta)
    % PARAMETERS
    s.N = length(F_tilde)-1;
    s.T = length(theta)-1;
    num_points = 100;
    F_vec = linspace(0, F_tilde(end), num_points);
    s.F_tilde = F_tilde;
    s.theta = theta;

    for k = 1:s.m
      for t = 1:s.T
        % FLUX CHARACTERISTIC 
        cvx_begin
          variables a(s.N) b(s.N)
          J = 0;
          for j = 1:s.N, for l = 1:num_points
            F = F_vec(l);
            if F_tilde(j) <= F && F <= F_tilde(j+1)
              J = J + square(s.f(F, theta(t), k) - a(j)*F - b(j));
            end
          end, end
          minimize(J);
          for j = 1:s.N-1
            a(j)*F_tilde(j+1) + b(j) == a(j+1)*F_tilde(j+1) + b(j+1);
          end
          b(1) == 0;
        cvx_end
        s.a{k}(:,t) = a;
        s.b{k}(:,t) = b;

        % PHASE TORQUE
        cvx_begin
          variables c(s.N) d(s.N)
          J = 0;
          for j = 1:s.N, for l = 1:num_points
            F = F_vec(l);
            if F_tilde(j) <= F && F <= F_tilde(j+1)
              J = J + square(s.g(F, theta(t), k) - c(j)*F - d(j));
            end
          end, end
          minimize(J);
          for j = 1:s.N-1
            c(j)*F_tilde(j+1) + d(j) == c(j+1)*F_tilde(j+1) + d(j+1);
          end
          d(1) == 0;
        cvx_end
        s.c{k}(:,t) = c;
        s.d{k}(:,t) = d;

      end
    end
  end

  function plot_pwa_3d(s)
    for k = 1:3, for t = 1:s.T
      f_tilde(t + (k-1)*s.T, 1) = s.b{k}(1,t);
      g_tilde(t + (k-1)*s.T, 1) = s.d{k}(1,t);
      for j = 1:s.N
        f_tilde(t + (k-1)*s.T, j+1) = s.a{k}(j,t)*s.F_tilde(j+1) + s.b{k}(j,t);
        g_tilde(t + (k-1)*s.T, j+1) = s.c{k}(j,t)*s.F_tilde(j+1) + s.d{k}(j,t);
      end
    end, end
    f_tilde = [f_tilde; f_tilde(1,:)];
    g_tilde = [g_tilde; g_tilde(1,:)];
    figure;
    theta = [s.theta(1:end-1), s.theta(1:end-1) + s.theta(end), s.theta + 2*s.theta(end)];
    mesh(theta, s.F_tilde, f_tilde')
    set(gca, 'XLim', [0, 2*pi])
    set(gca, 'YLim', [0, 6000]);
    set(gca, 'ZLim', [0, 4.5e-3]);
    print -depsc phi_hat.eps
    figure;
    mesh(theta, s.F_tilde, g_tilde')
    set(gca, 'XLim', [0, 2*pi])
    set(gca, 'YLim', [0, 6000]);
    set(gca, 'ZLim', [-80, 80]);
    print -depsc tau_hat.eps
  end


  function [I, v, tau, z, q, F, Psi, phi, lambda, time] = optimize(s, mu, tau_des, omega)
    T = s.T;
    l = s.l;
    m = s.m;
    n = s.n;
    N = s.N;
    dt = (s.theta(2) - s.theta(1))/omega;
    tic

    cvx_begin

      % VARIABLES
      cvx_solver gurobi
      variables F(m,T) I(n,T) z(m,T,N) v(n,T) Psi(m,T) phi(l,T) tau(m,T) lambda(n,T+1)
      variable q(m,T,N) binary

      % OBJECTIVE
      J = 0;
      if ~isequal(s.D, []); % use reformulated objective
        for t = 1:T, 
          for k = 1:m, for j = 1:N
            J = J + dt*s.D(k,k)*quad_over_lin(z(k,t,j), q(k,t,j));
          end, end
          J = J + dt*quad_form(F(:,t), s.Q - s.D);
        end 
      else % use standard objective
        for t = 1:T
          J = J + dt*quad_form(I(:,t), s.R);
        end
      end
      minimize(J + dt*sum(mu*square(sum(tau) - tau_des)))

      % CONSTRAINTS
      % voltage constraint
      abs(v) <= s.Vmax;

      % torque constraint
      mean(sum(tau(:,1:T))) == tau_des;

      % pwa variables
      sum(q, 3) == 1;
      sum(z, 3) == F;
      q >= 0;

      % periodicity
      lambda(3,1) == lambda(2,end);
      lambda(2,1) == lambda(1,end);
      lambda(1,1) == lambda(3,end);

      for t = 1:T

        % magnetic circuit
        s.M*F(:,t) == s.C*I(:,t);
        Psi(:,t) == s.M'*phi(:,t);
        lambda(:,t) == s.C'*phi(:,t);

        % dynamics
        lambda_prime = (lambda(:,t+1) - lambda(:,t))/(s.theta(t+1) - s.theta(t));
        v(:,t) == s.R*I(:,t) + omega*lambda_prime;

        % pwa approximations
        for k = 1:m
          sum_azbq = 0;
          sum_czdq = 0;
          for j = 1:N
            sum_azbq = sum_azbq + s.a{k}(j,t)*z(k,t,j) + s.b{k}(j,t)*q(k,t,j);
            sum_czdq = sum_czdq + s.c{k}(j,t)*z(k,t,j) + s.d{k}(j,t)*q(k,t,j);
            s.F_tilde(j)*q(k,t,j) <= z(k,t,j) <= s.F_tilde(j+1)*q(k,t,j);
          end
          Psi(k,t) == sum_azbq;
          tau(k,t) == sum_czdq;
        end
      end
    cvx_end
    time = toc;
  end

  function [theta, wave] = onephase2three(s, theta, wave)
    P = [0 1 0; 0 0 1; 1 0 0];
    wave = [ wave, P*wave, P^2*wave, wave(:,1)];
    theta = [ theta(1:end-1), theta(end) + theta(1:end-1), ... 
              2*theta(end) + theta(1:end-1), 3*theta(end)];
  end

  function s = find_D(s)
    s.Q = s.M'*(s.C\s.R)*(s.C\s.M);
    cvx_begin quiet
      variable D(s.m,s.m) diagonal
      maximize(trace(D))
      s.Q - D == semidefinite(s.m);
      D >= 0
    cvx_end
    s.D = D;
  end

end
end
