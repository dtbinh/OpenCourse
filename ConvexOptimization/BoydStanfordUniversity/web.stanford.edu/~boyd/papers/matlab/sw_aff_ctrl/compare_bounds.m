clear;
cvx_solver gurobi

% problem parameters
n = 2;
K = 3;
T = 10;
x_max = 5;
for i = 1:K
  A{i} = eye(n) + .1*randn(n);
  b{i} = .1*randn(n, 1);
  x0 = randn(n, 1);
  for j = 1:n
    M{i}(j,1) =  x_max * norm(A{i}(j,:), 1) + b{i}(j);
    m{i}(j,1) = -x_max * norm(A{i}(j,:), 1) + b{i}(j);
  end
end

% PERSPECTIVE-BASED FORMULATION
% global solution
[J_star_pers, x_pers, s_pers, runtime_pers] = sw_ctrl_pers(A, b, x0, T, x_max, 1);

% lower bound (integer-relaxation)
[lb_pers, ~, ~, ~] = sw_ctrl_pers(A, b, x0, T, x_max, 0);

% upper bound (relax-and-round)
[~, ~, s, ~] = sw_ctrl_pers(A, b, x0, T, x_max, 0);
[~, u] = max(s);
x_sim = x0;
for t_sim = 1:T
  x_sim(:,t_sim+1) = A{u(t_sim)}*x_sim(:,t_sim) + b{u(t_sim)};
end
ub_pers_rr = sum(sum(x_sim(:,2:end).^2));
if any(norms(x_sim(:,2:end), inf) > x_max);
  ub_pers_rr = inf;
end

% upper bound (shrinking-horizon heuristic)
x_sim = x0;
for t_sim = 1:T
  [~, ~, s, ~] = sw_ctrl_pers(A, b, x_sim(:,t_sim), T - t_sim + 1, x_max, 0);
  [~, u] = max(s(:, 1));
  x_sim(:,t_sim+1) = A{u}*x_sim(:,t_sim) + b{u};
end
ub_pers = sum(sum(x_sim(:,2:end).^2));
if any(norms(x_sim(:,2:end), inf) > x_max);
  ub_pers = inf;
end

% MIXED LOGICAL DYNAMICAL FORMULATION
% global solution
[J_star_mld, x_mld, s_mld, runtime_mld] = sw_ctrl_mld(A, b, x0, T, x_max, m, M, 1);

% lower bound (integer-relaxation)
[lb_mld, ~, ~, ~] = sw_ctrl_mld(A, b, x0, T, x_max, m, M, 0);

% upper bound (relax-and-round)
[~, ~, s, ~] = sw_ctrl_mld(A, b, x0, T, x_max, m, M, 0);
[~, u] = max(s);
x_sim = x0;
for t_sim = 1:T
  x_sim(:,t_sim+1) = A{u(t_sim)}*x_sim(:,t_sim) + b{u(t_sim)};
end
ub_mld_rr = sum(sum(x_sim(:,2:end).^2));
if any(norms(x_sim(:,2:end), inf) > x_max);
  ub_mld_rr = inf;
end

% upper bound (shrinking-horizon heuristic)
x_sim = x0;
for t_sim = 1:T
  [~, ~, s, ~] = sw_ctrl_mld(A, b, x_sim(:,t_sim), T - t_sim + 1, x_max, m, M, 0);
  [~, u] = max(s(:, 1));
  x_sim(:,t_sim+1) = A{u}*x_sim(:,t_sim) + b{u};
end
ub_mld = sum(sum(x_sim(:,2:end).^2));
if any(norms(x_sim(:,2:end), inf) > x_max);
  ub_mld = inf;
end


% GENERALIZED DISJUNCTIVE PROGRAMMING FORMULATION
% global solution
cvx_solver_settings('MSK_DPAR_MIO_MAX_TIME', 100000)
[J_star_gdp, x_gdp, s_gdp, runtime_gdp] = sw_ctrl_gdp(A, b, x0, T, x_max, 1);
cvx_solver_settings('MSK_DPAR_MIO_MAX_TIME', inf)

% lower bound (integer-relaxation)
[lb_gdp, ~, ~, ~] = sw_ctrl_gdp(A, b, x0, T, x_max, 0);

% upper bound (relax-and-round)
[~, ~, s, ~] = sw_ctrl_gdp(A, b, x0, T, x_max, 0);
[~, u] = max(s);
x_sim = x0;
for t_sim = 1:T
  x_sim(:,t_sim+1) = A{u(t_sim)}*x_sim(:,t_sim) + b{u(t_sim)};
end
ub_gdp_rr = sum(sum(x_sim(:,2:end).^2));
if any(norms(x_sim(:,2:end), inf) > x_max);
  ub_gdp_rr = inf;
end

% upper bound (shrinking-horizon heuristic)
x_sim = x0;
for t_sim = 1:T
  [~, ~, s, ~] = sw_ctrl_gdp(A, b, x_sim(:,t_sim), T - t_sim + 1, x_max, 0);
  [~, u] = max(s(:, 1));
  x_sim(:,t_sim+1) = A{u}*x_sim(:,t_sim) + b{u};
end
ub_gdp = sum(sum(x_sim(:,2:end).^2));
if any(norms(x_sim(:,2:end), inf) > x_max);
  ub_gdp = inf;
end

fprintf('\n\n')
fprintf('MLD relaxation lower bound is : %3.3f\n', lb_mld)
fprintf('MLD shrinking-horizon upper bound is : %3.3f\n', ub_mld)
fprintf('MLD relax-and-round upper bound is : %3.3f\n', ub_mld_rr)
fprintf('MLD formulation (global) solution : %3.3f\n', J_star_mld)
fprintf('\n')
fprintf('GDP relaxation lower bound is : %3.3f\n', lb_gdp)
fprintf('GDP shrinking-horizon upper bound is : %3.3f\n', ub_gdp)
fprintf('GDP relax-and-round upper bound is : %3.3f\n', ub_gdp_rr)
fprintf('GDP formulation (global) solution : %3.3f\n', J_star_gdp)
fprintf('\n')
fprintf('perspective relaxation lower bound is : %3.3f\n', lb_pers)
fprintf('perspective shrinking-horizon upper bound is : %3.3f\n', ub_pers)
fprintf('perspective relax-and-round upper bound is : %3.3f\n', ub_pers_rr)
fprintf('perspective formulation (global) solution : %3.3f\n', J_star_pers)
fprintf('\n')
