% Generates figure 11.7 of Boyd & Vandenberghe, Convex Optimization
%
% Barrier method for standard form LPs
% 
%     minimize   c'*x 
%     subject to A*x = b
%                x >= 0

function figure_11_7()

% Figure 11.7.  Gap versus Newton iterations for three problems.

disp('m=50, n=100.')
m = 50; 
n = 100;
A = randn(m,n);  
x0 = rand(n,1);  
b = A*x0; 
z = randn(m,1);  
c = A'*z + rand(n,1);
[x, iters, gaps] = stdlp(A,b,c,x0,100,1e-3,1e-5);
l = length(gaps);  iters1 = [];  gaps1 = [];
for i=1:l-1
    iters1 = [iters1,  iters(i)-1,  iters(i+1)-1];  
    gaps1 = [gaps1, gaps(i), gaps(i)]; 
end;
iters1 = [iters1, iters(l)-1] - iters1(1); 
gaps1 = [gaps1, gaps(l)]; 

disp('m=500, n=1000.')
m = 500; 
n = 1000;
A = randn(m,n); 
x0 = rand(n,1); 
b = A*x0; 
z = rand(m,1);  
c = A'*z + rand(n,1);
[x, iters, gaps] = stdlp(A,b,c,x0,100,1e-3,1e-5);
l = length(gaps);  iters2 = [];  gaps2 = [];
for i=1:l-1
   iters2 = [iters2, iters(i)-1, iters(i+1)-1];  
   gaps2 = [gaps2, gaps(i), gaps(i)]; 
end;
iters2 = [iters2, iters(l)-1] - iters2(1);  
gaps2 = [gaps2, gaps(l)]; 

disp('m=1000, n=2000.')
m = 1000; 
n = 2000;
A = randn(m,n); 
x0 = rand(n,1); 
b = A*x0; 
z = rand(m,1);  
c = A'*z + rand(n,1);
[x, iters, gaps] = stdlp(A,b,c,x0,100,1e-3,1e-5);
l = length(gaps);  iters3 = [];  gaps3 = [];
for i=1:l-1
   iters3 = [iters3, iters(i)-1, iters(i+1)-1];  
   gaps3 = [gaps3, gaps(i), gaps(i)]; 
end;
iters3 = [iters3, iters(l)-1] - iters3(1);  
gaps3 = [gaps3, gaps(l)]; 


semilogy(iters1,gaps1, iters2,gaps2, '--', iters3,gaps3, '-');
text(iters1(length(iters1)), gaps1(length(iters1)),'m=50');
text(iters2(length(iters2)), gaps2(length(iters2)),'m=100');
text(iters3(length(iters3)), gaps3(length(iters3)),'m=500');
xlabel('Newton iterations'); ylabel('Duality gap');


function [x, inniters, gaps] = stdlp(A, b, c, x0, mu, tol, reltol)

% [x, inniters, gaps] = stdlp(A, b, c, x0, mu, tol, reltol)
%
% (primal) minimize  c'*x        (dual)  maximize   b'*z 
%          subjec to A*x = b             subject to A'*z <= c
%                    x >= 0
%
% x0:  strictly feasible starting point, not necesarily on central path
% inniters:  array with number of Newton iterations per outer iteration
% gaps:  array with duality gaps at the end of each outer iteration

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
    val = t*c'*x - sum(log(x));
    g = t*c - 1./x;
    z =  -(A*spdiags(x.^2,0,n,n)*A') \ (A*(g.*(x.^2)) );
    dx = - (g+A'*z) .* (x.^2);
    fprime = (g+A'*z)'*dx;
    if ((-fprime/2) < NTTOL)
        z = -z/t;
        dualobj = b'*z;
        gap = x'*(c-A'*z);  
        inniters = [inniters, k];  
        gaps = [gaps,gap];
        if (gap < tol),  return;  end;
        t = min(t*mu, (n+1)/tol);  
    else
        s = 1;
        while (min(x+s*dx) < 0), s = BETA*s; end;
        while (s*(t*c+A'*z)'*dx - sum(log(1+s*dx./x)) > ...
            ALPHA*s*fprime), 
            s = BETA*s; 
        end;
        x = x+s*dx;
    end;
end;
