% Generates figures 11.15 and 11.16 of Boyd & Vandenberghe, Convex 
% Optimization
%
% Barrier method  for a small SOCP.

function figure_11_15()

% Generate a pair of strictly primal and dual feasible SOCPs.
%
% (primal) minimize    c'*x
%          subject to  ||Ai*x + bi||_2 <= ci'*x + di, i=1,...,m
%
% (dual)   maximize    b'*z 
%          subject to  A'*z = c
%                      ||zi||_2 <= wi, i=1,...,m
%
% where 
%
% A = [A1; c1'; A2; c2'; ..., Am; cm']
% b = [b1; d1; b2; d2; ..., bm; dm]
% z = [z1; w1; z2; w2; ... ; zm; wm]

randn('state',0);
rand('state',0);
m = 50; 
n = 50;  
szs = 5*ones(m,1);   % row dimensions of the Ai's
N = sum(szs);

% indsi are the indices of ci' in A
% indsf are the indices of the first row of Ai in A 
% E is a matrix such that [sum(y1); ... ; sum(ym) ] = E*y
indsi = cumsum(szs+1);
indsf = cumsum(szs+1)-szs;
E = sparse(m,N+m);
for i=1:m
    E(i,indsf(i):indsi(i)-1) = ones(1,szs(i));
end;

% Ransom x0;  random zi's;  wi = ||zi|| + random positive slack.
x0 = randn(n,1);  
z = randn(N+m,1);  
z(indsi) = sqrt(E*(z.^2)) + rand(m,1);
A = randn(N+m,n);   
c = A'*z;
b = randn(N+m,1);
b(indsi) = 0;
y = A*x0+b;
b(indsi) = sqrt(E*(y.^2)) - y(indsi) + rand(m,1);

% Compute optimal value and shift b0 so that optimal value is 1.
[x, inniters, gaps] = socp(c, A, b, szs, x0, 1e-6, 10);
s = (1-c'*x)/norm(c)^2;
b = b - s*A*c;
x0 = x0+s*c; 


% xc is the point on the central path with barrier parameter t=1.
% To center, call socp() with tol = 2*m+1.

[xc, inniters, gaps] = socp(c, A, b, szs, x0, 2*m+1, 10);

% Figure 11.15.   Gap versus Newton iterastions for three values of mu.

figure(1)
muvals = [2, 50, 200];
[x, iters, gaps] = socp(c, A, b, szs, xc, 1e-6, muvals(1));
l = length(gaps);  iters1 = [];  gaps1 = [];
for i=1:l-1
    iters1 = [iters1, iters(i)-1, iters(i+1)-1];  
    gaps1 = [gaps1, gaps(i), gaps(i)]; 
end;
iters1 = [iters1, iters(l)-1];  gaps1 = [gaps1, gaps(l)]; 

[x, iters, gaps] = socp(c, A, b, szs, xc, 1e-6, muvals(2));
l = length(gaps);  iters2 = [];  gaps2 = [];
for i=1:l-1
    iters2 = [iters2, iters(i)-1, iters(i+1)-1];  
    gaps2 = [gaps2, gaps(i), gaps(i)]; 
end;
iters2 = [iters2, iters(l)-1];  gaps2 = [gaps2, gaps(l)]; 

[x,iters,gaps] = socp(c, A, b, szs, xc, 1e-6, muvals(3));
l = length(gaps);  iters3 = [];  gaps3 = [];
for i=1:l-1
    iters3 = [iters3, iters(i)-1, iters(i+1)-1];  
    gaps3 = [gaps3, gaps(i), gaps(i)]; 
end;
iters3 = [iters3, iters(l)-1];  gaps3 = [gaps3, gaps(l)]; 

semilogy(iters1,gaps1, iters2,gaps2, '-', iters3,gaps3, '--');
text(iters1(length(iters1)), gaps1(length(iters1)),'mu=2');
text(iters2(length(iters2)), gaps2(length(iters2)),'mu=50');
text(iters3(length(iters3)), gaps3(length(iters3)),'mu=200');
axis([0 90  0.5e-7 1e3]);
xlabel('Newton iterations');  ylabel('duality gap');


% Figure 11.16.  Number of Newton iterations versus mu.

muvals = [1.5, 2:2:8, 10:10:200];
noiters=zeros(1,length(muvals));
for i=1:length(muvals)
    [x, iters, gaps] = socp(c, A, b, szs, xc, 1e-3, muvals(i));
    noiters(i) = iters(length(iters));
end;
figure(2)
plot(muvals,noiters,'o', muvals,noiters,'-');
axis([ 0 200 0 150]);
xlabel('mu');  ylabel('Newton iterations');


function [x, inniters, gaps] = socp(c, A, b, szs, x0, tol, mu)

%
% [x, inniters, gaps] = socp(c, A, b, szs, x0, tol, mu)
%
% Solves the SOCP
%
%     minimize    c'*x
%     subject to  ||Ai*x+bi||^2/(ci'*x+di) <= ci'*x+di, i=1,...,m
%                 ci'*x + di >= 0                           
%
% given a strictly primal feasible x0.
%
% A = [A1; c1'; A2'; c2'; ... ; Am; cm']
% b = [b1; d1; b2; d2; ...; bm; dm]
% szs = [n1; n2; ... ; nm] where Ai = nixm
%
% The dual problem is
%
%     maximize    b'*z 
%     subject to  A'*z = c
%                 ||zi|| <= wi
%              
%  with variable z = [z1; w1; z2; w2; ... ; zm; wm].


MAXITERS = 200;  
NTTOL = 1e-5;    % stop centering when Nt decrement <= NTTOL 
ALPHA = 0.01; 
BETA = 0.5;      % parameters for backtracking line search

n = length(x0);  
m = length(szs);   
N = sum(szs);

% indsi are the indices of ci' in A
% indsf are the indices of the first row of Ai in A 
% E is a matrix s.t. [sum(y1);  ... ; sum(ym) ]' = E*y
% E2 is a matrix s.t. [ sum(y1)-w1;  ... ' sum(ym)-wm ]' = E2*y

indsi = cumsum(szs+1);  
indsf = cumsum(szs+1)-szs;
E = sparse(m,N+m);  
E2 = sparse(m,N+m);
for i=1:m
    E(i,indsf(i):indsi(i)-1) = ones(1,szs(i));
    E2(i,indsf(i):indsi(i)) = [ones(1,szs(i)) -1];
end;

x = x0;  
y = A*x+b;  
primalobj = c'*x;  
dualobj = -Inf;  
gap = Inf;
t = 1;  
inniters = [];  gaps = []; 
for iters = 1:MAXITERS
    y = A*x+b;  
    v = y(indsi);  
    s = v.^2 - E*(y.^2);

    % phi = t*c'*x + sum_i f(ci'*x+di, Ai*x+bi) where  
    % f(w,v) = - sum_i log(v^2 - w'*w) 
    % grad f(w,v) = (2/(v^2 - w'*w)) * [w; -v];
    % hess f(w,v) = (2*(v^2-w'*w)*I + 4*[w;-v]*[w;-v]')/(v^2 - w'*w)^2
    val = t*c'*x - sum(log(s));
    d = 2*(y./(E2'*s));  
    g = t*c + A'*d;
    H = A' * (spdiags(2./(E2'*s),0,N+m,N+m) * A);
    for i=1:m
        u = A(indsf(i):indsi(i),:)'*d(indsf(i):indsi(i));
        H = H + u*u';
    end;
    dx = -H\g;
    fprime = g'*dx;
    dy = A*dx;
    s = 1;
    newx = x+s*dx;  newy = A*newx+b;  newv = newy(indsi);
    while (min(newv < 0) | (max(sqrt(E*(newy.^2)) - newv) > 0)) 
        s = BETA*s; 
        newx = x+s*dx;  newy = A*newx+b;  newv = newy(indsi);
    end;
    newval = t*c'*newx - sum(log(newv)) - ...
        sum(log(newv- (E*newy.^2) ./ newv));
    while (newval > val + s*ALPHA*fprime)
        s = BETA*s; 
        newx = x+s*dx;  newy = A*newx+b;  newv = newy(indsi);
        newval = t*c'*newx - sum(log(newv)) - ...
            sum(log(newv- (E*newy.^2) ./ newv));
    end;
    x = x+s*dx;
    if (-fprime/2 < NTTOL), 
        gap = 2*m/t; 
        gaps = [gaps, gap];
        inniters = [inniters, iters]; 
        if (gap < tol), return; end;
        t = min(t*mu, (2*m+1)/tol);  
    end; 
end;
