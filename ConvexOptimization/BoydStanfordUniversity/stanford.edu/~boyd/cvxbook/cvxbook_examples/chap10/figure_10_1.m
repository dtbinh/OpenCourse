% Generates figures 10.1, 10.2, 10.3, 10.4 of Boyd & Vandenberghe,
% Convex Optimization
%
% Infeasible start Newton method for solving
%
%     minimize   -sum_i log x_i
%     subject to Ax = b


% Generate a feasible problem.

rand('state',3);   
randn('state',4);
n = 100;        
m = 50;        
A = randn(m,n);
x0 = rand(n,1);
b = A*x0;        


% Apply infeasible Newton method to KKT conditions
%
%     -1./x + A'*w = 0 
%      A*x = b

alpha = 0.011;
beta  = 0.5;
tol = 1e-12;
maxiters = 20;

x = ones(n,1);  
w = zeros(m,1); 
resp = zeros(1, maxiters);
resd = zeros(1, maxiters);
steps = zeros(1, maxiters);
for k=1:maxiters
    rp = A*x-b;  rd = -1./x + A'*w;  
    resp(k) = norm(rp);  resd(k) = norm(rd);
    normr = norm([rd; rp]);
    if normr < 1e-12, break; end;
    dw = - (A * diag(x.^2) * A') \ (A*((x.^2).*rd) - rp); 
    dx = - (x.^2) .* (rd + A'*dw);
    s = 1;
    while min(x+s*dx) <= 0,  s = s*beta; end;
    while norm([ -1./(x+s*dx) + A'*(w+s*dw);  A*(x+s*dx)-b ]) > ...
        (1-alpha*s) * normr
        s = s*beta;
    end;
    steps(k) = s;
    x = x + s*dx;
    w = w + s*dw;
end;

% Figure 10.1
figure(1)
semilogy([0:k-1], resp(1:k), 'o', [0:k-1], resp(1:k), '-', ...
    [0:k-1], resd(1:k), 'o', [0:k-1], resd(1:k), '--');
xlabel('iteration number')
ylabel('residuals')

% Figure 10.2 
figure(2);
plot([1:k-1], steps(1:k-1), 'o', [1:k-1], steps(1:k-1), '-');
xlabel('iteration number')
ylabel('step length')
axis([0 12 0.4 1.1]);


% Generate an infeasible problem.

rand('state',3);
randn('state',4);
A = randn(m,n);
x0 = randn(n,1); 
b = A*x0;       


% Apply infeasible Newton method to KKT conditions
%
%     -1./x + A'*w = 0 
%      A*x = b

alpha = 0.011;
beta  = 0.5;
tol = 1e-12;
maxiters = 21;

x = ones(n,1);  
w = zeros(m,1); 
resp = zeros(1, maxiters);
resd = zeros(1, maxiters);
steps = zeros(1, maxiters);
for k=1:maxiters
    rp = A*x-b;  rd = -1./x + A'*w;  
    resp(k) = norm(rp);  resd(k) = norm(rd);
    normr = norm([rd; rp]);
    if normr < 1e-12, break; end;
    dw = - (A * diag(x.^2) * A') \ (A*((x.^2).*rd) - rp); 
    dx = - (x.^2) .* (rd + A'*dw);
    s = 1;
    while min(x+s*dx) <= 0,  s = s*beta; end;
    while norm([ -1./(x+s*dx) + A'*(w+s*dw);  A*(x+s*dx)-b ]) > ...
        (1-alpha*s) * normr
        s = s*beta;
    end;
    steps(k) = s;
    x = x + s*dx;
    w = w + s*dw;
end;

% Figure 10.3
figure(3)
semilogy([0:k-1], resp(1:k), 'o', [0:k-1], resp(1:k), '-', ...
    [0:k-1], resd(1:k), 'o', [0:k-1], resd(1:k), '--');
xlabel('iteration number')
ylabel('residuals')

% Figure 10.4 
figure(4);
plot([1:k-1], steps(1:k-1), 'o', [1:k-1], steps(1:k-1), '-');
xlabel('iteration number')
ylabel('step length')
axis([0 20 0 0.3]);
