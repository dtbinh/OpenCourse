% Generates figures 11.4 and 11.5 of Boyd & Vandenberghe, Convex 
% Optimization
%
% Barrier method  for a small LP.

function figure_11_4()

% Generate a strictly primal and dual feasible LP.
%
% (primal) minimize    c'*x 
%          subject to  A*x <= b
% 
% (dual)   maximize    -b'*z 
%          subject to  A'*z + c = 0, z >= 0

randn('state',0);
rand('state',0);
m = 100;
n = 50;
A = randn(m,n);
b = rand(m,1);
c = A'*rand(m,1); 

% Solve using linprog.
x = linprog(c, A, b);
optval = c'*x;

% Make a change of variables  x := x + s*c, so that optimal value is 1.
s = (1-optval)/(c'*c);
b = b + s*A*c;
x0 = s*c; 


% xc is the point on the central path with barrier parameter t=1.
t = 1;  
x = x0;
for k=1:100
    d = b-A*x;
    val = t*c'*x - sum(log(d));
    g = t*c + A'*(1./d);
    H = A'*diag(1./(d.^2))*A;
    v = -H\g;  
    fprime = g'*v;
    s = 1;
    while (min(b-A*(x+s*v)) < 0),  s = s/2;  end;
    while (t*c'*(x+s*v) - sum(log(b-A*(x+s*v))) > val + 0.01*s*fprime), 
        s = s/2; 
    end;
    x = x+s*v;
    if ((-fprime/2) < 1e-6), break; end;
end;
xc = x;


% Figure 11.4.  Barrier method starting at xc, for three values of mu.

muvals = [2 50 150];
[x,iters,gaps] = lp(A,b,c,xc,muvals(1),1e-6);
l = length(gaps);  iters1 = [];  gaps1 = [];
for i=1:l-1
    iters1 = [iters1 iters(i)-1 iters(i+1)-1];  
    gaps1 = [gaps1 gaps(i) gaps(i)]; 
end;
iters1 = [iters1 iters(l)-1]; gaps1 = [gaps1 gaps(l)]; 

[x,iters,gaps] = lp(A,b,c,xc,muvals(2),1e-6);
l = length(gaps);  iters2 = [];  gaps2 = [];
for i=1:l-1
    iters2 = [iters2 iters(i)-1 iters(i+1)-1];  
    gaps2 = [gaps2 gaps(i) gaps(i)]; 
end;
iters2 = [iters2 iters(l)-1]; gaps2 = [gaps2 gaps(l)]; 

[x,iters,gaps] = lp(A,b,c,xc,muvals(3),1e-6);
l = length(gaps);  iters3 = [];  gaps3 = [];
for i=1:l-1
    iters3 = [iters3 iters(i)-1 iters(i+1)-1];  
    gaps3 = [gaps3 gaps(i) gaps(i)]; 
end;
iters3 = [iters3 iters(l)-1]; gaps3 = [gaps3 gaps(l)]; 

figure(1)
semilogy(iters1,gaps1, iters2,gaps2, '-', iters3,gaps3, '--');
axis([0 90 0.5e-7 1e3]);
text(iters1(length(iters1)), gaps1(length(iters1)),'mu=2');
text(iters2(length(iters2)), gaps2(length(iters2)),'mu=50');
text(iters3(length(iters3)), gaps3(length(iters3)),'mu=100');
axis([0 90 0.5e-7 1e3]);
xlabel('Newton iterations');  ylabel('duality gap');


% Figure 11.5.  Number of Newton iterations vs. mu.

muvals = [1.5 2 4 6 8 10:10:200];
noiters=zeros(1,length(muvals));
for i=1:length(muvals)
    [x,iters,gaps] = lp(A,b,c,xc,muvals(i),1e-3);
    noiters(i) = iters(length(iters));
end;
figure(2)
plot(muvals,noiters,'o', muvals,noiters,'-');
axis([ 0 200 0 150]);
xlabel('mu'); ylabel('Newton iterations');



function [x, inniters, gaps] = lp(A, b, c, x0, mu, tol)

% [x, inniters, gaps] = lp(A, b, c, x0, mu, tol)
%
%      minimize    c'*x 
%      subject to  A*x <= b
%
%      maximize    -b'*z
%      subject to  A'*z + c = 0
%                  z >= 0
%
% Barrier method for solving LP within absolute accuracy tol, starting
% with initial t = 1, at a strictly feasible x0.  We assume the 
% problem is strictly dual feasible.
% 
% inniters:  array with number of Newton iters per outer iteration
% gaps:   array with duality gaps at the end of each outer iteration
 

MAXITERS = 500;   
ALPHA = 0.01;
BETA = 0.5;
NTTOL = 1e-5;     % stop inner iteration if lambda^2/2 < NTTOL

[m,n] = size(A);
t = 1;
x = x0;
gaps = [];
inniters = [];
for k=1:MAXITERS
    d = b-A*x;
    val = t*c'*x - sum(log(d));
    g = t*c + A'*(1./d);
    H = A'*diag(1./(d.^2))*A;
    v = -H\g;  
    fprime = g'*v;
    s = 1;
    while (min(b-A*(x+s*v)) < 0),  s = BETA*s;  end;
    while (t*c'*(x+s*v) - sum(log(b-A*(x+s*v))) > val + ALPHA*s*fprime),
        s = BETA*s; 
    end;
    x = x+s*v;
    if ((-fprime/2) < NTTOL)
        inniters = [inniters, k];  
        z = (1./d) .* (1 + (A*v)./d) / t;
        gap = (b-A*x)'*z;  
        gaps = [gaps,gap];
        if (gap < tol), return;  end; 
        t = min(t*mu, (m+1)/tol);  
   end;
end;
inniters = [inniters, k];  
gaps = [gaps, gap];  
disp(['Maxiters (', int2str(MAXITERS), ') exceeded.']);
