% Generates figures 9.7 and 9.8 of Boyd & Vandenberghe, Convex
% Optimization
%
% Effect of scaling on convergence of gradient method for the problem
%
%     minimize c'*x - sum_i log(bi-ai'*x)  
%
% x is scaled as T*x with 
%
%     T = diag(1, gamma^(1/n), gamma^(2/n), ..., gamma^((n-1)/n)

randn('state',0);
rand('state',0);   
A = randn(500,100);
b = rand(500,1);         
c = 20*randn(100,1);
m = 500;
n = 100;


% Compute minimum using Newton's method.

alpha = 0.01;
beta = 0.5;
maxiters = 50;
x = zeros(100,1);
for k=1:maxiters
    y = b-A*x;
    val = c'*x - sum(log(y));  
    grad = c + A'*(1./y);
    hess = A'*diag(1./(y.^2))*A;
    if (norm(grad) < 1e-8), break; end;
    v = -hess\grad; 
    fprime = grad'*v;
    s = 1;
    dy = -(A*v)./y;
    dc = c'*v;
    for lsiters = 1:100
        newx = x+s*v;
        if (min(1+s*dy) < -1e-5) 
            s = beta*s;
        else 
            newval = val + s*dc - sum(log(1+s*dy));
            if (newval > val + s*alpha*fprime), s = beta*s;
            else break;
            end;
        end;
    end;
    x =  x+s*v;
end;
fmin = val;
xmin = x;


% Figure 9.8.  Condition numbers at optimum vs gamma.

figure(1)
gammas = logspace(-1.4,1.5,1000);
condnrs = zeros(size(gammas));
for k=1:length(gammas);
   gamma = gammas(k);
   T = diag(gamma.^[0:1/n:(n-1)/n]);
   condnrs(k) = cond(T'*hess*T);
end;
loglog(gammas, condnrs, '-');
axis([0.05 20 10 1e4]); 
xlabel('gamma');
ylabel('kappa');


% Figure 9.7.  Number of iterations versus gamma.

gammas  = logspace(-1.5,1.5,11);
iters = zeros(size(gammas));
alpha = 0.3;
beta = 0.7;
maxiters = 5000;
for i = 1:length(gammas)
    gamma = gammas(i);
    T = diag(gamma.^[0:1/n:(n-1)/n]);
    x = zeros(n,1);
    Asc = A*T;
    csc = T'*c;
    for k=1:maxiters
        y = b-Asc*x;
        val = csc'*x - sum(log(y));  
        grad = csc + Asc'*(1./y);
        if (val - fmin < 1e-5), break; end;
        v = -grad; 
        fprime = grad'*v;
        s = 1;
        dy = -(Asc*v)./y;
        dc = csc'*v;
        for lsiters = 1:1000
            if (min(1+s*dy) < -1e-5) 
                s = beta*s;
            else 
                newval = val + s*dc - sum(log(1+s*dy));
                if (newval > val + s*alpha*fprime), s = beta*s;
                else break;
                end;
            end;
        end;
    x =  x+s*v;
    end;
    iters(i) = k;
end;

figure(2)
loglog(gammas, iters, 'o', gammas, iters,'-');
axis([0.05 20 100 5000]); 
xlabel('gamma');
ylabel('iters');
