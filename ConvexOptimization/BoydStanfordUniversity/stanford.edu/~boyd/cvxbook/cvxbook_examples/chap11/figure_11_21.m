% Generates figure 11.21 of Boyd & Vandenberghe, Convex Optimization
%
% Primal-dual method for a small LP.

function figure_11_21()

% Generate a problem with optimal value 1.

randn('state',0);
rand('state',0);
m = 100;
n = 50;
A = randn(m,n);
b = rand(m,1);
z0 = rand(m,1);
c = -A'*z0;
x = linprog(c, A, b);
s = (1-c'*x)/norm(c)^2;
b = b + s*A*c;
x0 = s*c; 

[x, z, iters, gaps, rs] = lp_pd(A, b, c, x0);

figure(1)
semilogy([0:length(gaps)-1], gaps,'-', [0:length(gaps)-1], gaps,'o');
xlabel('iteration number');  ylabel('surrogate gap');

figure(2)
semilogy([0:length(gaps)-1], rs, '-', [0:length(gaps)-1], rs, 'o');
axis([0 30 1e-16 1e5]);
xlabel('iteration number'); ylabel('residual');


function [x, z, iters, gaps, resdls] = lp_pd(A, b, c, x0)

% [x, z, iters, gaps, resdls] = lp_pd(A, b, c, x0)
%
% Solves 
%
%     minimize    c'*x                 maximize    -b'*z
%     subject to  A*x <= b             subject to  A'*z + c = 0
%                                                  z >= 0 
% 
% using the primal-dual method with simple line search on norm
% of residual.
%
% x0 must be primal feasible.


MAXITERS = 500;  
TOL = 1e-8;       
RESTOL = 1e-8;    
MU = 10;
ALPHA = 0.01;
BETA = 0.5;

[m,n] = size(A);
gaps = [];  resdls = [];
x = x0;  s = b-A*x;  z = 1./s;  
for iters = 1:MAXITERS
    gap = s'*z;  gaps = [gaps, gap];
    res = A'*z  + c;   resdls = [resdls, norm(res)];
    if ((gap < TOL) & (norm(res) < RESTOL)), return;  end;
    tinv = gap/(m*MU);
    sol = -[zeros(n,n), A'; A, diag(-s./z)] \ [A'*z+c; -s+tinv*(1./z)];
    dx = sol(1:n);  dz = sol(n+[1:m]);  ds = -A*dx;

    % backtracking line search
    r = [c+A'*z;  z.*s-tinv];
    step = min(1.0, 0.99/max(-dz./z));
    while (min(s+step*ds) <= 0), step = BETA*step; end;
    newz = z+step*dz;  newx = x+step*dx;  news = s+step*ds;
    newr = [c+A'*newz;  newz.*news-tinv];
    while (norm(newr) > (1-ALPHA*step)*norm(r))
        step = BETA*step;
        newz = z+step*dz;  newx = x+step*dx;  news = s+step*ds;
        newr = [c+A'*newz;  newz.*news-tinv];
    end;
    x = x+step*dx;  z = z +step*dz;  s = b-A*x;
end;
