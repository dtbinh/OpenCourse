% Frobenious norm diagonal scaling example
% Section 4.5.4 example in Boyd & Vandenberghe "Convex Optimization"
% (see page 163 for more details)
%
% Given a square matrix M, the goal is to find a positive vector d
% such that ||DMD^{-1}||_F is minimized, where D = diag(d).
% The problem can be cast into an unconstrained GP:
%
%   minimize   sum_{i,j=1}^{n} M_ij^2*d_i^2/d_j^2
%
% with variable d (vector in R^n) and data matrix M in R^{n-by-n}.
%
% Almir Mutapcic 10/15/05
randn('state',0);

% matrix size (M is an n-by-n matrix)
n = 4;
M = randn(n,n);

% GP variables
gpvar d(n);

% objective expressed using a for loop
obj = posynomial; % initialize an empty posynomial
for i = 1:n
  for j = 1:n
    obj = obj + M(i,j)^2*d(i)^2*d(j)^-2;
  end
end

% objective can also be expressed using matrices
% (can also use it below in the problem solving routine)
% obj_vect = sum( sum( diag(d.^2)*(M.^2)*diag(d.^-2) ) );

% solve the unconstrained GP problem
[opt_frob_norm, solution, status] = gpsolve(obj,[]);

% assign GP variables to their values (d now becomes an array of doubles)
assign(solution);

% construct matrix D and display results
disp(' ')
disp('Optimal diagonal scaling d is: '), d
D = diag(d);

frob_norm_M = norm(M,'fro');
frob_norm_A = norm(D*M*inv(D),'fro');
fprintf(1,['The Frobenius norm for M before the scaling is %3.4f\n'...
          'and after the optimal scaling (DMD^{-1}) it is %3.4f.\n'],...
          frob_norm_M,frob_norm_A);
