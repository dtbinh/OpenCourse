% Generates figure 9.23 of Boyd & Vandenberghe, Convex Optimization
%
% Newton method for a problem with 10000 variables.

% minimize c'*x - sum_i log(bi-ai'*x);  
randn('state',0);
rand('state',0);   
m = 100000;  
n = 10000;
N = 200;
A = [ speye(n,n); -speye(n,n);  sprandn(m-2*n,n,0.000001) ];
m1 = (m-2*n)/N;  n1 = n/N;    
A1 = sprandn(m1,n1,0.05);
A = A + [sparse(2*n,n); kron(speye(N,N),sparse(A1))];
b = [ones(n,1); ones(n,1); rand(m-2*n,1)];  
figure(1); 
spy(A'*A);  
title('sparsity pattern of hessian');

maxiters = 50;
alpha = 0.01;
beta = 0.5;
vals = zeros(1,maxiters);
steps= zeros(1,maxiters);
x = zeros(n,1);
for iter=1:maxiters
    y = b-A*x;  yinv = 1./y;
    val = -sum(log(y));  vals(iter) = val;
    grad = A'*yinv;
    D = spdiags(yinv.^2, 0, m, m);  
    hess = A'*D*A;
    hess1 = (hess+hess')/2;
    v = -hess1\grad; 
    lambdasq = -grad'*v;
    disp([int2str(iter), ': Nt. decr. = ', num2str(sqrt(lambdasq))]);
    if (sqrt(lambdasq) < 1e-5), break; end;
    fprime = -lambdasq;
    s = 1;
    dy = -(A*v)./y;
    while (min(1+s*dy) < 0), s = beta*s; end;
    newx = x+s*v;  newval = val - sum(log(1+s*dy));
    while (newval > val + s*alpha*fprime)
        s = beta*s;
        newx = x+s*v;  
        newval = val - sum(log(1+s*dy));
    end;
    x =  x+s*v;
end;
fmin = vals(iter);

figure(2)
semilogy([0:iter-1], vals(1:iter) - fmin, '-', ...
         [0:iter-1], vals(1:iter) - fmin, 'o');
axis([0 20 1e-8 1e5]);
xlabel('k')
ylabel('error')
