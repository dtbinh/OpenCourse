% CVX examples for ICM 2006 talk available at
% http://www.stanford.edu/~boyd/cvx_opt_graph_lapl_eigs.html
%
% This is an example with a small graph (8 nodes and 13 edges).
% Written for CVX by Almir Mutapcic 08/29/06

% small example (incidence matrix A)
A = [ 1  0  0  1  0  0  0  0  0  0  0  0  0;
     -1  1  0  0  1  1  0  0  0  0  0  0  1;
      0 -1  1  0  0  0  0  0 -1  0  0  0  0;
      0  0 -1  0  0 -1  0  0  0 -1  0  0  0;
      0  0  0 -1  0  0 -1  1  0  0  0  0  0;
      0  0  0  0  0  0  1  0  0  0  1  0  0;
      0  0  0  0  0  0  0 -1  1  0 -1  1 -1;
      0  0  0  0 -1  0  0  0  0  1  0 -1  0];

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

%% look at the weighted Laplacian eigenvalues
%eig(eye(n) - A*diag(w_fdla)*A' - 1/n*ones(n))
%eig(eye(n) - A*diag(w_fmmc)*A' - 1/n*ones(n))
%eig(eye(n) - A*diag(w_md)*A'   - 1/n*ones(n))
%eig(eye(n) - A*diag(w_bc)*A'   - 1/n*ones(n))
%eig(eye(n) - A*diag(w_mh)*A'   - 1/n*ones(n))

fprintf(1,'\nResults:\n');
fprintf(1,'FDLA weights:\t\t rho = %5.4f \t tau = %5.4f\n',rho_fdla,tau_fdla);
fprintf(1,'FMMC weights:\t\t rho = %5.4f \t tau = %5.4f\n',rho_fmmc,tau_fmmc);
fprintf(1,'M-H weights:\t\t rho = %5.4f \t tau = %5.4f\n',rho_mh,tau_mh);
fprintf(1,'MAX_DEG weights:\t rho = %5.4f \t tau = %5.4f\n',rho_max_deg,tau_max_deg);
fprintf(1,'BEST_CONST weights:\t rho = %5.4f \t tau = %5.4f\n', ...
        rho_best_const,tau_best_const);
