% Generates figure 11.6 of Boyd & Vandenberghe, Convex Optimization.
% 
% Barrier method for a small GP.

function figure_11_6()

% Generate a pair of strictly primal and dual feasible GPs
%
% (primal)  minimize    f(A0*x + b0)
%           subject to  f(Ai*x + bi) <= 0, i = 1,...,m
%
% (dual)    maximize    b0'*nu0 + ... + bm'*num + g(nu0, 1) + 
%                       g(nu1, sum(nu1)) + ...  + g(num,sum(num))
%           subject to  sum(nu0) = 1
%                       A0'*nu0 + ... + Am'*num = 0
%                       nui >= 0, i = 0,...,m
%
% where f(y) = log(sum(exp(y),  g(y,z) = -sum(y.*log(y/z))

rand('state',0);
randn('state',0);
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

[x, iters, gaps] = gp(A, b, szs, x0, 10, 1e-6);
optval = log(sum(exp(A(1:m0,:)*x + b(1:m0))));
b(1:m0) = b(1:m0) + 1 - optval;

% xc is the point on the central path with barrier parameter t=1.
% To center with t=1, call gp() with tol = m+1.

[xc, iters, gaps] = gp(A, b, szs, x0, 20, m+1);


% Figure 11.6.  Gap versus Newton iterations for three values of mu.

figure(1)
muvals = [2, 50, 150];
[x, iters, gaps] = gp(A, b, szs, xc, muvals(1), 1e-6);
l = length(gaps);  iters1 = [];  gaps1 = [];
for i=1:l-1
    iters1 = [iters1, iters(i)-1, iters(i+1)-1];  
    gaps1 = [gaps1, gaps(i), gaps(i)]; 
end;
iters1 = [iters1, iters(l)-1];  gaps1 = [gaps1, gaps(l)]; 

[x, iters, gaps] = gp(A, b, szs, xc, muvals(2), 1e-6);
l = length(gaps);  iters2 = [];  gaps2 = [];
for i=1:l-1
    iters2 = [iters2, iters(i)-1, iters(i+1)-1];  
    gaps2 = [gaps2, gaps(i), gaps(i)]; 
end;
iters2 = [iters2, iters(l)-1];  gaps2 = [gaps2, gaps(l)]; 

[x, iters, gaps] = gp(A, b, szs, xc, muvals(3), 1e-6);
l = length(gaps);  iters3 = [];  gaps3 = [];
for i=1:l-1
    iters3 = [iters3, iters(i)-1, iters(i+1)-1];  
    gaps3 = [gaps3, gaps(i), gaps(i)]; 
end;
iters3 = [iters3, iters(l)-1];  gaps3 = [gaps3, gaps(l)]; 

semilogy(iters1,gaps1, iters2,gaps2, '-', iters3,gaps3, '--');
text(iters1(length(iters1)), gaps1(length(iters1)),'mu=2');
text(iters2(length(iters2)), gaps2(length(iters2)),'mu=50');
text(iters3(length(iters3)), gaps3(length(iters3)),'mu=150');
axis([0 120  0.5e-7 1e3]);
xlabel('Newton iterations');  ylabel('duality gap');


function [x, inniters, gaps] = gp(A, b, szs, x0, mu, abstol)

% [x, inniters, gaps] = gp(A, b, szs, x0, mu, abstol)
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
% given a strictly primal feasible x0.
%
% f is defined as f(y) = log(sum(exp(y)) 
% g is defined as g(y,z) = -sum(y.*log(y/z))
%
% Input arguments.
% - A:        (sum_i n_i)xn matrix; A = [A0' A1' ... Am']'
% - b:        (sum_i n_i) vector; b = [b0' b1' ... bm' ]'
% - szs:      dimensions of Ai and bi;  szs = [n0 n1 ... nm]' where Ai 
%             is (nixn) and bi is (nix1)
% - x0:       n-vector; must be strictly feasible
% - abstol:   absolute tolerance 
%
% Output arguments.
% - x:         n-vector; primal optimum. 
% - inniters:  array with Newton iterations per outer iteration.
% - gaps:      array with gaps after each outer iteration.


MAXITERS = 500;  
NTTOL = 1e-5;     % stop centering when Nt decrement <= NTTOL 
ALPHA = 0.01;     % backtracking parameters
BETA = 0.5;      

n = length(x0);  
N = size(A,1);  
m = length(szs)-1;
n0 = szs(1);

% indsf is an array with indices of the first rows of each block in A
% indsl is an array with indices of the last rows of each block in A
% E is a matrix s.t. [sum(y0), sum(y1),  ..., sum(ym) ]' = E*y
indsf = cumsum(szs)-szs+1;  
indsl= cumsum(szs);         
E = sparse(m+1,N);
for i=1:m+1
    E(i,indsf(i):indsl(i)) = ones(1,szs(i));
end;

x = x0;  
y = A*x+b;  
f = log(E*exp(y));
t = 1;  

ntiters = 0;
inniters = [];
gaps = [];
while (ntiters <= MAXITERS)


    % Centering: minimize phi(x) 
    %
    % phi(x)  = t*f(A0*x+b0) - sum_i log(-f(Ai*x+bi))  
    % gradphi = t*A0'*gradf(y0)) + sum_i (-1/f(yi)) Ai'*gradf(yi)
    % hessphi = t*A0'*hessf(y0)*A0 + ...
    %           sum_i Ai' * ((1/f(yi))^2 *gradf(yi)*gradf(yi)' + ...
    %           (-1/f(yi))*hessf(yi)) * Ai
    %
    % where
    %
    % gradf(y) = exp(y) / sum(exp(y))
    % hessf(y) = -(exp(y)*exp(y)') / sum(exp(y))^2 + ...
    %            diag(exp(y))/sum(exp(y))

    while (ntiters <= MAXITERS)
        gradf = exp(y)./(E'*(E*exp(y)));
        gradphi = A' * (gradf .* (E'*[t;-1./f(2:m+1)]));
        hessphi = ...   % first the gradf terms
            (A' * (sparse(diag(gradf)) * E' * ...
            sparse(diag([0;1./f(2:m+1).^2])))) * ...
            ((E * sparse(diag(gradf))) * A)  ...
            - ... % add the rank one hessf terms
            (A' * (sparse(diag(exp(y))) * E' * ...
            sparse(diag([t; -1./f(2:m+1)]./(E*exp(y)).^2)))) * ...
            ((E * sparse(diag(exp(y)))) * A)  ... 
            + ... % add the diagonal hessf terms
            (A'*sparse(diag(gradf.*(E'* [t;-1./f(2:m+1)])))) * A;
        dx = -hessphi\gradphi;
        ntiters = ntiters + 1;
        fprime = gradphi'*dx;
        dy = A*dx;
        if (-fprime/2 < NTTOL), break; end; 
        phi = t*f(1) - sum(log(-f(2:m+1)));
        s = 1;
        for lsiters = 1:50
            newf = log(E*exp(y+s*dy));
            if (max(newf(2:m+1)) < 0)
                newphi = t*newf(1) - sum(log(-newf(2:m+1)));
                if (newphi < phi + ALPHA*s*fprime), break; end;
            end;
            s = BETA*s;
        end;
        x = x+s*dx;
        y = A*x+b;
        f = log(E*exp(y));
    end;

    % Dual variables on central path, corrected to account for 
    % incomplete centering.
    %
    % nu0 = (1 - g0'*dy0)*g0 + diag(dy0)*g0 
    % nui = (-1/t*fi) * ( (1-(1+1/fi)*gi'*dyi)*gi + diag(dyi)*gi )
    % lambdai = (-1/tfi) * (1-gi'*dyi/fi) 

    lambda = -(1-E(2:m+1,:)*(gradf.*dy)./f(2:m+1)) ./ (t*f(2:m+1));
    nu =  [ (1-gradf(1:n0)'*dy(1:n0) + dy(1:n0)) .* gradf(1:n0);
        ( E(2:m+1,n0+1:N)' * ( 1- ((1+1./f(2:m+1)) .* (E(2:m+1,:) * ...
        (gradf.*dy))) ) + dy(n0+1:N) ) .* gradf(n0+1:N) ./ ...
        (-t * E(2:m+1, n0+1:N)'*f(2:m+1)) ]; 
    dualobj = b'*nu - sum(nu.*log(nu./(E'*E*nu)));
    primalobj = f(1);
    gap = primalobj-dualobj; 
    inniters = [inniters, ntiters];
    gaps = [gaps, gap];
    if (gap < abstol), return; end;
    t = min(mu*t, (m+1)/abstol);  

end;
