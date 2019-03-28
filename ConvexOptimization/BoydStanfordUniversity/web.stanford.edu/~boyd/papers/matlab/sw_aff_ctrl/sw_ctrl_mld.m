function [J_star, x, s, cpu_time] = sw_ctrl_mld(A, b, x0, T, x_max, m, M, binary);
n = length(x0);
K = length(A);

tic

cvx_begin
  
  % decision variables
  variables y(n, T, K) x(n, T+1)
  if binary, variable s(K, T) binary; else, variable s(K, T); end

  % objective
  J = 0;
  for t = 1:T
    J = J + sum_square(x(:,t+1));
  end
  minimize(J);

  % dynamics and logic
  for t = 1:T, for i = 1:K,
    if t == 1
      y(:,1,i) == s(i,1)*(A{i}*x0 + b{i});
    else
      A{i}*x(:,t) + b{i} - M{i}*(1-s(i,t)) <= y(:,t,i);
      y(:,t,i) <= A{i}*x(:,t) + b{i} - m{i}*(1-s(i,t));
      m{i}*s(i,t) <= y(:,t,i) <= M{i}*s(i,t);
    end
  end, end
  x(:,2:T+1) == sum(y, 3);
  sum(s, 1) == 1;
  0 <= s <= 1;

  % state constraints
  -x_max <= x(:,2:T+1) <= x_max

  % initial condition
  x(:,1) == x0;

cvx_end

cpu_time = toc;
J_star = cvx_optval;
