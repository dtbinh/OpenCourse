% Generates figure 10.5 of Boyd & Vandenberghe, Convex Optimization
%
% Newton's method applied to a convex-concave game
%
%     f(x,y) = x'*A*y + c'*x + d'*y - log(1-x'*x) + log(1-y'*y);


randn('state',0);
n = 100;
A = randn(n,n);
c = randn(n,1);
d = randn(n,1);

BETA = .5;
ALPHA = .01;
x = .01*ones(n,1);
y = .01*ones(n,1);
resx = []; 
resy = []; 
res = []; 

% Infeasible start Newton method for solving r(x,y) = 0 where 
%
%     r(x,y) = ( A*y + c + 2*x/(1-x'*x),  A'*x + d - 2*y/(1-y'*y) )   
 
for iters =1:100

    r = [ A*y + (2/(1-x'*x))*x + c;  A'*x - (2/(1-y'*y))*y + d ];   
    resx = [resx, norm(r(1:n))];
    resy = [resy, norm(r(n+[1:n]))];
    res = [res, norm(r)];
    if (norm(r) < 1e-8), break;  end; 

    Dr = [ ((2/(1-x'*x))*eye(n) + (4/(1-x'*x)^2)*x*x'), A ;
        A',  (-(2/(1-y'*y))*eye(n) - (4/(1-y'*y)^2)*y*y') ];
    step = -Dr\r;  dx = step(1:n);  dy = step(n+[1:n]);

    s = 1;
    newx = x+s*dx;  newy = y+s*dy;
    while ((norm(newx) >= 1) | (norm(newy) >= 1)), 
        s = BETA*s;  
        newx = x+s*dx;   newy = y+s*dy;
    end;
    newr = [ A*newy + (2/(1-newx'*newx))*newx + c; 
        A'*newx - (2/(1-newy'*newy))*newy + d ];  
    while (norm(newr) > (1-ALPHA*s)*norm(r))
        s = BETA*s;  
        newx = x+s*dx;   newy = y+s*dy;
        newr = [ A*newy + (2/(1-newx'*newx))*newx + c; 
            A'*newx - (2/(1-newy'*newy))*newy + d];   
    end;
    x = x+s*dx;  y = y+s*dy;

end;

semilogy(0:(length(res)-1), res, 'o', 0:(length(res)-1), res, '-');
xlabel('iteration');
ylabel('residual');
axis([0 8 1e-15 1e5]);
