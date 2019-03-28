function [J_star, x, s, cpu_time] = sw_ctrl_pers(A, b, x0, T, x_max, binary);
n = length(x0);
K = length(A);

tic

cvx_begin %quiet
  
  % decision variables
  variables z(n, T, K) x(n, T+1)
  if binary, variable s(K, T) binary; else, variable s(K, T); end

  % objective
  J = 0;
  for t = 1:T
    J = J + sum_square(x(:,t+1));
  end
  minimize(J);

  % dynamics
  for t = 0:T-1
    Azb_sum = 0;
    for i = 1:K
      Azb_sum = Azb_sum + A{i}*z(:,t+1,i) + b{i}*s(i,t+1);
    end
    x(:,t+2) == Azb_sum;
  end

  % other constraints
  0 <= s <= 1;
  sum(s, 1) == 1;
  x(:,1:T) == sum(z, 3);

  % state constraint
  for t = 1:T-1, for i = 1:K
    abs(z(:,t+1,i)) <= s(i,t+1)*x_max;
  end, end
  abs(x(:,T+1)) <= x_max;

  % initial condition
  for i = 1:K
    z(:,1,i) == x0*s(i,1);
  end

cvx_end

cpu_time = toc;
J_star = cvx_optval;
