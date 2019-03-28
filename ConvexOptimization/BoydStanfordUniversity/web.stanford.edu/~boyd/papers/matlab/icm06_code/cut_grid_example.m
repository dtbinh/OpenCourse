% CVX examples for ICM 2006 talk available at
% http://www.stanford.edu/~boyd/cvx_opt_graph_lapl_eigs.html
%
% This is an example with a cut-grid graph (64 nodes and 95 edges).
% Written for CVX by Almir Mutapcic 08/29/06

%********************************************************************
% generate a cut-grid graph
%********************************************************************
% cut grid example (incidence matrix A)
[A,B,xy] = cut_grid_data;
A = sparse(A);

%********************************************************************
% compute edge weights
%********************************************************************
fprintf(1,'WARNING: The optimal weight computations take some time...\n');
% size of the network
[n,m] = size(A);

% compute various edge weights (some optimal, some based on heuristics)
w_fdla = fdla(A);
w_fmmc = fmmc(A);
w_md   = max_deg(A);
w_bc   = best_const(A);
w_mh   = mh(A);

% report results
rho_fdla       = norm(eye(n) - A*diag(w_fdla)*A'  - 1/n*ones(n));
rho_fmmc       = norm(eye(n) - A*diag(w_fmmc)*A'  - 1/n*ones(n));
rho_max_deg    = norm(eye(n) - A*diag(w_md)*A'  - 1/n*ones(n));
rho_best_const = norm(eye(n) - A*diag(w_bc)*A'  - 1/n*ones(n));
rho_mh         = norm(eye(n) - A*diag(w_mh)*A'  - 1/n*ones(n));

tau_fdla       = 1/log(1/rho_fdla);
tau_fmmc       = 1/log(1/rho_fmmc);
tau_max_deg    = 1/log(1/rho_max_deg);
tau_best_const = 1/log(1/rho_best_const);
tau_mh         = 1/log(1/rho_mh);

eig_fdla       = eig(eye(n) - A*diag(w_fdla)*A' - 1/n*ones(n));
eig_fmmc       = eig(eye(n) - A*diag(w_fmmc)*A' - 1/n*ones(n));
eig_max_deg    = eig(eye(n) - A*diag(w_md)*A'   - 1/n*ones(n));
eig_best_const = eig(eye(n) - A*diag(w_bc)*A'   - 1/n*ones(n));
eig_mh         = eig(eye(n) - A*diag(w_mh)*A'   - 1/n*ones(n));

fprintf(1,'\nResults:\n');
fprintf(1,'FDLA weights:\t\t rho = %5.4f \t tau = %5.4f\n',rho_fdla,tau_fdla);
fprintf(1,'FMMC weights:\t\t rho = %5.4f \t tau = %5.4f\n',rho_fmmc,tau_fmmc);
fprintf(1,'M-H weights:\t\t rho = %5.4f \t tau = %5.4f\n',rho_mh,tau_mh);
fprintf(1,'MAX_DEG weights:\t rho = %5.4f \t tau = %5.4f\n',rho_max_deg,tau_max_deg);
fprintf(1,'BEST_CONST weights:\t rho = %5.4f \t tau = %5.4f\n', ...
        rho_best_const,tau_best_const);

%********************************************************************
% plot results
%********************************************************************
figure(1), clf
plotgraph(A,xy,w_fdla);

figure(2), clf
plotgraph(A,xy,w_fmmc);
