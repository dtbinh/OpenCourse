% Generates figures 9.6, 9.21, 9.22 of Boyd & Vandenberghe, Convex
% Optimization
%
% Gradient descent and Newton's method for a problem with 100 variables.

% f(x) = c'*x - sum_i log(bi-ai'*x);  
randn('state',0);
rand('state',0);   
A = randn(500,100);
b = rand(500,1);         
c = 20*randn(100,1);

% Compute minimum using Newton's method.

alpha = 0.01 ;
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
            else, break;
            end;
        end;
    end;
    x =  x+s*v;
end;
fmin = val;

% Figure 9.6.  Gradient method with backtracking line search.

figure(1)
alpha = 0.1;
beta = 0.5;
maxiters = 2000;
x = zeros(100,1);
vals = zeros(1,maxiters);
steps= zeros(1,maxiters);
for k=1:maxiters
    y = b-A*x;
    val = c'*x - sum(log(y));  vals(k) = val;
    grad = c + A'*(1./y);
    if (val-fmin < 1e-4), break; end;
    v = -grad; 
    fprime = grad'*v;
    s = 1;
    dy = -(A*v)./y;
    dc = c'*v;
    for lsiters = 1:1000
        if (min(1+s*dy) < -1e-5), 
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
semilogy([0:k-1], vals(1:k) - fmin, '-');
text(150, vals(150)-fmin,'backtracking line search');
xlabel('k');
ylabel('error');
hold on

% Figure 9.6: Gradient method with exact line search.

figure(1)
x = zeros(100,1);
vals = zeros(1,maxiters);
steps= zeros(1,maxiters);
for k=1:maxiters
    y = b-A*x;
    val = c'*x - sum(log(y));  vals(k) = val;
    grad = c + A'*(1./y);
    if (val-fmin < 1e-4), break; end;
    v = -grad; 
    % exact line search by bisection
    dy = -(A*v)./y;
    dc = c'*v;
    l = 0;   gl = v'*grad;  
    inds = find(dy<0);
    u = 0.99*min(-1./dy(inds));
    gr = dc - sum(dy./(1+u*dy));
    ustart = u;
    while ((u-l)/ustart > 1e-3)
        s = (u+l)/2;
        gs = dc - sum(dy./(1+s*dy));
        if (gs<0), l=s; 
        elseif (gs>0), u=s; 
        else, break, end;
    end;
    x =  x+s*v;
end;
semilogy([0:k-1], vals(1:k) - fmin, '--');
text(50,vals(50)-fmin,'exact line search');
xlabel('k');
ylabel('error');
axis([0 200 1e-4 1e4]);
hold off


% Figures 9.21 and 9.22.  Newton's method with backtracking.

alpha = 0.01 ;
beta = 0.5;
maxiters = 50;
x = zeros(100,1);
vals = zeros(1,maxiters);
steps= zeros(1,maxiters);
for k=1:maxiters
    y = b-A*x;
    val = c'*x - sum(log(y));  vals(k) = val;
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
            else, break;
            end;
        end;
    end;
    steps(k) = s;
    x =  x+s*v;
end;

% Figure 9.21.  Error vs. iteration.
figure(2)
semilogy([0:k-1], vals(1:k) - fmin, '-', [0:k-1], vals(1:k) - fmin, ...
    'o');
text(k-3,vals(k-2)-fmin,'backtracking line search');
hold on


% Figure 9.22.  Step size vs. iteration.
figure(3)
plot([1:k-1], steps(1:k-1), '-', [1:k-1], steps(1:k-1), 'o');
text(2,steps(2),'backtracking line search');
hold on


% Figures 9.21 and 9.22.  Newton's method with exact line search.

maxiters = 50;
x = zeros(100,1);
vals = zeros(1,maxiters);
steps= zeros(1,maxiters);
for k=1:maxiters
    y = b-A*x;
    val = c'*x - sum(log(y));  vals(k) = val;
    grad = c + A'*(1./y);
    hess = A'*diag(1./(y.^2))*A;
    if (norm(grad) < 1e-8), break; end;
    v = -hess\grad; 
    dy = -(A*v)./y;
    dc = c'*v;
    l = 0;   gl = v'*grad;  
    inds = find(dy<0);
    u = 0.99*min(-1./dy(inds));
    gr = dc - sum(dy./(1+u*dy));
    ustart = u;
    while ((u-l) > 1e-3)
        s = (u+l)/2;
        gs = dc - sum(dy./(1+s*dy));
        if (gs<0), l=s; 
        elseif (gs>0), u=s; 
        else, break, end;
    end;
    steps(k) = s;
    x =  x+s*v;
end;

% Error vs. iteration.  Add to figure  9.21.
figure(2)
semilogy([0:k-2], vals(1:k-1) - fmin, '-', [0:k-2], ...
    vals(1:k-1) - fmin, 'd');
text(k-3,vals(k-2)-fmin,'exact line search');
hold off
xlabel('k')
ylabel('error')
axis([0 10 1e-15 1e5]);
xlabel('k'); ylabel('error')
hold off

% Step vs. iteration.  Add to figure 9.22.
figure(3)
plot([1:k-1], steps(1:k-1), '-', [1:k-1], steps(1:k-1), 'd');
text(3,steps(3),'exact line search');
axis([0 8 0 2.0]);
hold off
