% Generates figures 10.6, 10.7, 10.8 of Boyd & Vandenberghe, Convex
% Optimization
%
% Three methods for equality constrained analytic centering
%
%     (primal) minimize    -sum(log(x))    
%              subject to  A*x = b
%
%     (dual)   maximize   -b'*z + sum(log(A'*z)) + n 


% Generate random primal and dual feasible problem

randn('state',3);
rand('state',3);
m = 100;  n = 500;
A = randn(m,n);

% make z0 dual feasible (A'*z0 = 1)
z0 = randn(m,1); 
A = A + z0*(1 - A'*z0)'/(z0'*z0);   

% make x0 primal feasible
x0 = rand(n,1);   
b = A*x0;   


MAXITERS = 50;
ALPHA = 0.1;
BETA  = 0.5;
TOL = 1e-12;


% Figure 10.6.   Feasible Newton method applied to primal problem, 
% with four different starting points.

figure(1); 
randn('state',2);
nostarts = 4;

for starts = 1:nostarts

    % Starting point: generate random vector v in nullspace of A and 
    % take x = x0 + 0.99*s*v where s is maximum step to the boundary.
    v = randn(n,1);   v = v - A'*(A'\v);  
    x = x0 + (0.99/max(-v./x0))*v;

    fvals = [];
    for iters=1:MAXITERS
        f = -sum(log(x));  fvals = [fvals, f];
        g = -1./x;  Hinv = x.^2;
        dw = - ( (A.*(Hinv(:,ones(1,m)))')*A') \  (A*(Hinv.*g)); 
        dx = - Hinv.* (g+A'*dw);
        lambdasq = -g'*dx;
        if (lambdasq/2 < TOL), break; end;
        s = 1; 
        while (min(x+s*dx) < 0), s = BETA*s; end;
        while (-sum(log(1+x+s*dx)) > f - ALPHA*lambdasq), 
            s = BETA*S; 
        end;
        x = x+s*dx;
    end;

    figure(1);
    semilogy([0:iters-2], fvals(1:(iters-1))-fvals(iters), 'o', ...
        [0:iters-2], fvals(1:(iters-1))-fvals(iters), 'b-');
    hold on

end;

axis([0 20 1e-10 1e5]);
xlabel('iteration');
ylabel('error');
hold off


% Figure 10.7.  Newton's method applied to the dual.
% 
%     mininimize b'*z - sum(log(A'*z)

figure(2); 
randn('state',1);
nostarts = 4;

for starts = 1:nostarts

    % Starting point: generate random vector v and take 
    % z = z0 + 0.99*s*v where s is maximum step to the boundary. 
    v = randn(m,1);   
    z = z0 + (0.99/max(-(A'*v)./(A'*z0)))*v;

    fvals = [];
    for iters=1:MAXITERS
        y = A'*z;
        f = b'*z - sum(log(y));  fvals = [fvals, f];
        g = b + -A*(1./y); 
        d = 1./(y.^2);
        H = (A .* (d(:,ones(1,m))'))*A';
        dz = -H\g;
        lambdasq = -g'*dz;
        if (lambdasq/2 < TOL), break; end;
        dy = A'*dz;
        s = 1; 
        while (min(y+s*dy) < 0), s = BETA*s; end;
        while (b'*(z+s*dz)-sum(log(A'*(z+s*dz))) > f-ALPHA*lambdasq) 
	    s = BETA*s; 
        end;
        z = z+s*dz;
    end;

    figure(2);
    semilogy([0:iters-2], fvals(1:(iters-1)) - fvals(iters), 'o', ...
        [0:iters-2], fvals(1:(iters-1)) - fvals(iters), 'b-');
        hold on;

end;

xlabel('iteration');  ylabel('error');
axis([0 10 1e-10 1e5]);
hold off


% Figure 1.8.  Infeasible Newton  method.

figure(3); 
randn('state',2);
nostarts = 4;

for starts = 1:nostarts

    x = rand(n,1);  
    w = randn(m,1);
    resdls = [];   
    for iters=1:MAXITERS
        g = -1./x;  
        r = [g + A'*w;  A*x-b];  
        normr = norm(r);  resdls = [resdls; normr];
        if (normr) < sqrt(TOL), break; end;
        H = 1./(x.^2);  Hinv = x.^2;
        v = ((A.*(Hinv(:,ones(1,m)))')*A') \  (2*A*x-b);
        dw = v-w;
        dx = x - x.^2 .* (A'*v);
        s = 1; 
        while (min(x+s*dx) < 0), s = BETA*s; end;
        while norm([-1./(x+s*dx) + A'*(w+s*dw);  A*(x+s*dx)-b ]) > ...
            (1-ALPHA*s)*normr, 
            s = BETA*s; 
        end;
        x = x+s*dx;  w = w+s*dw;
    end;

   semilogy([0:iters-1], resdls(1:iters), 'o', [0:iters-1], ...
       resdls(1:iters), 'b-');
   hold on

end;

axis([0 25 1e-15 1e10]);
xlabel('iteration');
ylabel('residual');
axis([0 25 1e-15 1e10]);
hold off
