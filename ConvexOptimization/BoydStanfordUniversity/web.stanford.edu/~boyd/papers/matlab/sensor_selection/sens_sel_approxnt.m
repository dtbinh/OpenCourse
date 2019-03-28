function [zhat L ztilde Utilde] = sens_sel_approxnt(A, k)
% Solves the problem
%	maximize log det (sum_{i=1}^m z_i a_i a_i^T) + kappa sum_{i=1}^m(log(z_i)+ log(1-z_i))
%	subject to sum(z) = k
%			   0 <= z_i <= 1, i=1,..., m
% variable z in R^m
% problem parameters kappa (>0), a_1, ..., a_m in R^n
%
% see paper Sensor Selection via Convex Optimization
% www.stanford.edu/~boyd/papers/sensor_selection.html
%
% Nov 2007 Siddharth Joshi & Stephen Boyd

% Newton's method parameters
MAXITER  = 30;
NT_TOL = 1e-3;
GAP = 1.005;
% Backtracking line search parameters
alpha = 0.01;
beta = 0.5;

[m n] = size(A);
z = ones(m,1)*(k/m); % initialize
g = zeros(m,1);
ones_m = ones(m,1);
kappa = log(GAP)*n/m; 
% guarantees GM of lengths of semi-axes of ellipsoid corresponding to 
% ztilde <= 1.01 from optimal

fprintf('\nIter.  Step_size  Newton_decr.  Objective  log_det\n');

fz = -log(det(A'*diag(z)*A)) - kappa*sum(log(z) + log(1-z));

fprintf('   0\t  -- \t     --   %10.3f  %10.3f\n', -fz, log(det(A'*diag(z)*A)));

for iter=1:MAXITER

    W = inv(A'*diag(z)*A);
    V = A*W*A';

    g = -diag(V)- kappa*(1./z - 1./(1-z));
    H = V.*V + kappa*diag(1./(z.^2) + 1./((1-z).^2));

    R = chol(H);
    Hinvg = (R\(R'\g));
    Hinv1 = (R\(R'\ones_m));
    dz = -Hinvg + ((ones_m'*Hinvg) / (ones_m'*Hinv1))*Hinv1;

    deczi = find(dz < 0);
    inczi = find(dz > 0);
    s = min([1; 0.99*[-z(deczi)./dz(deczi) ; (1-z(inczi))./dz(inczi)]]);

    while (1)
        zp = z + s*dz;
        fzp = -log(det(A'*diag(zp)*A)) - kappa*sum(log(zp) + log(1-zp));

        if (fzp <= fz + alpha*s*g'*dz)
            break;
        end
        s = beta*s;
    end
    z = zp; fz = fzp;
    
    fprintf('%4d %10.3f %10.3f %10.3f %10.3f\n', iter, s, -g'*dz/2, -fz, log(det(A'*diag(z)*A)));

    if(-g'*dz/2 <= NT_TOL)
        break;
    end
end

zsort=sort(z); thres=zsort(m-k); zhat=(z>thres);
L = log(det(A'*diag(zhat)*A));
ztilde = z; 
Utilde = log(det(A'*diag(z)*A)) + 2*m*kappa;
