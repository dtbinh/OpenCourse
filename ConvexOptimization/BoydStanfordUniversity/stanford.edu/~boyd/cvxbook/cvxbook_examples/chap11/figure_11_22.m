% Generates figure 11.22 of Boyd & Vandenberghe, Convex Optimization
%
% Primal-dual method for a small GP.

function figure_11_22()

% Generate a problem with optimal value 1

randn('state',0);
rand('state',0);
m = 100; 
n = 50;
szs = 5*ones(m+1,1);
N = sum(szs);
m0 = szs(1);

% Random postive nu0 with sum(nu0) = 1.
nu = rand(N,1);
nu(1:szs(1)) = nu(1:szs(1)) / sum(nu(1:szs(1)));

% Random A that satisfies A'*nu = 0.
A = randn(N,n);
A = A - nu*(nu'*A)/(nu'*nu);

% Random x0.
x0 = randn(n,1);

% Random b0.   Other bi's are constructed so that f(Ai*x0+bi) + zi = 0
% with random positive initial slacks z.
b = zeros(N,1);
b(1:m0) = randn(m0,1);
z = rand(m,1);
ind = m0;
for i = 1:m
    mi = szs(i+1);
    b(ind+(1:mi)) = -A(ind+(1:mi),:)*x0 - z(i) - log(mi);
    ind = ind + mi;
end

% Compute optimal value and shift b0 so that optimal value is 1.

[x, z, iters, gaps, resdls] = gp_pd(A, b, szs, x0);
optval = log(sum(exp(A(1:m0,:)*x + b(1:m0))));
b(1:m0) = b(1:m0) + 1 - optval;


[x, z, iters, gaps, rs] = gp_pd(A, b, szs, x0);

figure(1)
semilogy([0:length(gaps)-1], gaps,'-', [0:length(gaps)-1], gaps,'o');
xlabel('iteration number');  ylabel('surrogate gap');

figure(2)
semilogy([0:length(resdls)-1], rs, '-', [0:length(rs)-1], rs, 'o');
axis([0 length(resdls) 1e-15 1e5]);
xlabel('iteration number');  ylabel('residual');


function [x, z, iters, gaps, resdls] = gp_pd(A, b, szs, x0)

% [x, z, iters, gaps, resdls] = gp_pd(A, b, szs, x0)
%
% Solves the geometric program
%
%     minimize    f(y0)
%     subject to  f(yi) <= 0, i=1,...,m
%                 Ai*x+bi = yi, i=0,...,m
%
% and its dual
%
%     maximize    b0'*nu0 + ... + bm'*num + g(nu0, 1) + ...
%                 g(nu1,lambda1) + ,..  + g(num, lambdam)
%     subject to  sum(nu0) = 1
%                 sum(nui) = lambdai, i=1,...,m
%                 A0'*nu0 + ... + Am'*num = 0
%                 nui >= 0, i=0,...,m
%                 lambdai >= 0, i=1,...,m
%
% where f(y) = log(sum(exp(y),  g(y,z) = -sum(y.*log(y/z))
% using a primal-dual method with a simplified line search.

% Input arguments:
% - A:      (sum_i n_i)xn matrix; A = [A0; A1; ...; Am]
% - b:      (sum_i n_i) vector; b = [b0; b1; ...; bm]
% - szs:    array with row dimensions of Ai and bi
% - x0:     n-vector; must be strictly feasible

MAXITERS = 500;  
TOL = 1e-8;       
RESTOL = 1e-8;    
MU = 10;
ALPHA = 0.01;
BETA = 0.5;

n = length(x0);  
N = size(A,1);  
m = length(szs)-1;  
n0 = szs(1);

% indsf is an array with indices of the first rows of each block in A
% indsl is an array with indices of the last rows of each block in A
% E is a matrix s.t. [sum(y0), sum(y1),  ..., sum(ym) ]' = E*y
E = sparse(m+1,N);
indsf = cumsum(szs)-szs+1;  
indsl= cumsum(szs);        
for i=1:m+1
    E(i,indsf(i):indsl(i)) = ones(1,szs(i));
end;

x = x0;  y = A*x+b;  
f = log(E*exp(y));  s = -f(2:m+1);  optval = f(1);
z = 1./s;   
gaps = [];  resdls = [];
for iters = 1:MAXITERS
    y = A*x+b;  f = log(E*exp(y));  s = -f(2:m+1);  optval = f(1);
    gap = s'*z;
    res = A' * ( exp(y) .* (E'* ([1;z] ./ (E*exp(y)))) );
    w = exp(y) .* (E' * (1./(E*exp(y))));   
    gradsf = A' * (diag(w)*E');    % [gradf0 gradf1 ... gradfm]
    gaps = [gaps gap];
    resdls = [resdls norm(res)];
    if ((gap < TOL) & (norm(res) < RESTOL)), return;  end;
    H = A(indsf(1):indsl(1),:)' * ( diag(w(indsf(1):indsl(1))) - ...
       w(indsf(1):indsl(1)) * w(indsf(1):indsl(1))') * ...  
       A(indsf(1):indsl(1),:);
    for i=1:m
        H  = H + z(i) * A(indsf(i+1):indsl(i+1),:)' * ...
            ( diag(w(indsf(i+1):indsl(i+1))) - ...
            w(indsf(i+1):indsl(i+1)) * w(indsf(i+1):indsl(i+1))') * ... 
            A(indsf(i+1):indsl(i+1),:);
    end;
    tinv = gap/(m*MU);
    sol = -[ H, gradsf(:,2:m+1);  gradsf(:,2:m+1)' diag(-s./z) ] \ ...
        [ res ;  -s + tinv*(1./z) ];
    dx = sol(1:n);  dz = sol(n+[1:m]);  dy = A*dx;

    r = [res; z.*s-tinv];
    step = min(1.0, 0.99/max(-dz./z));
    newx = x + step*dx;  newy = y + step*dy;
    news = -log(E(2:m+1,:)*exp(newy));
    while (min(news) <= 0.0)
        step = BETA*step;
        newx = x + step*dx;  newy = y + step*dy;
        news = -log(E(2:m+1,:)*exp(newy));
    end;
    newz = z + step*dz;
    newres = A' * ( exp(newy) .* (E'* ([1;newz] ./ (E*exp(newy)))) );
    newr = [newres; newz.*news-tinv];
    while (norm(newr) > (1-ALPHA*step)*norm(r))
        step = BETA*step;
        newx = x + step*dx;  newy = y + step*dy;
        newz = z + step*dz;  news = -log(E(2:m+1,:)*exp(newy));
        newres = A' * (exp(newy) .* (E'* ([1;newz] ./ (E*exp(newy)))));
        newr = [newres; newz.*news-tinv];
    end;
    x = x+step*dx;  z = z +step*dz;
end;
