% Generates figure 11.19 of Boyd & Vandenberghe, Convex optimization
%
% Barrier method for SDP
%
%     minimize    sum(x)
%     subject to  B + diag(x) >= 0


function figure_11_19()

% Figure 11.19.  Gap versus Newton iterations for three problems.

randn('state',0);
rand('state',0);
mu = 20;

disp('n=50.')
n = 50;  
B = tril(randn(n,n));  
B = B + B' - 2*diag(diag(B));  
B = B/max(abs(eig(B)));
[x, z, iters, gaps] = etp(B, mu, 1e-4);
l = length(gaps);  iters1 = [];  gaps1 = [];
for i=1:l-1
    iters1 = [iters1, iters(i)-1, iters(i+1)-1];  
    gaps1 = [gaps1, gaps(i), gaps(i)]; 
end;
iters1 = [iters1, iters(l)-1] - iters1(1);  
gaps1 = [gaps1, gaps(l)]; 

disp('n=500.')
n = 500;  
B = tril(randn(n,n));  
B = B + B' - 2*diag(diag(B));  
B = B/max(abs(eig(B)));
[x, z, iters, gaps] = etp(B, mu, 1e-4);
l = length(gaps);  iters2 = [];  gaps2 = [];
for i=1:l-1
    iters2 = [iters2, iters(i)-1, iters(i+1)-1];  
    gaps2 = [gaps2, gaps(i), gaps(i)]; 
end;
iters2 = [iters2, iters(l)-1] - iters2(1);  
gaps2 = [gaps2, gaps(l)]; 

disp('n=1000.')
n = 1000;  
B = tril(randn(n,n));  
B = B + B' - 2*diag(diag(B)); 
B = B/max(abs(eig(B)));
[x, z, iters, gaps] = etp(B, mu, 1e-4);
l = length(gaps);  iters3 = [];  gaps3 = [];
for i=1:l-1
    iters3 = [iters3, iters(i)-1, iters(i+1)-1];  
    gaps3 = [gaps3, gaps(i), gaps(i)]; 
end;
iters3 = [iters3, iters(l)-1] - iters3(1);  
gaps3 = [gaps3, gaps(l)]; 

semilogy(iters1, gaps1, iters2, gaps2, '-', iters3, gaps3, '-');
xlabel('Newton iterations');  ylabel('duality gap');
text(iters1(length(iters1)), gaps1(length(iters1)),'n=50');
text(iters2(length(iters2)), gaps2(length(iters2)),'n=500');
text(iters3(length(iters3)), gaps3(length(iters3)),'n=1000');
axis([0 50 1e-5 1e5]);


function [x, Z, iters, gaps] = etp(B,mu,tol)

% [x, Z, iters, gaps] = etp(B, mu, tol)
%
% Solves the SDP
%
%     minimize    sum(x)
%     subject to  B + diag(x) >= 0
%
% B:  nxn symmetric
% iters:  array with number of Newton iterations per outer iteration
% gaps:  array with duality gaps at the end of each outer iteration
%
% The dual problem is
%
%     maximize   -tr(B*Z)
%     subject to diag(Z) = 1 
%                Z = 0
%
%

n = size(B,1);  

MAXITERS = 500;            
ALPHA = 0.01;  
BETA = 0.5;  
NTTOL = 1e-5;               % stop Nt if lambda^2/2 < NTTOL

gaps = [];  iters = [];
x = 1.1*(-min(eig(B)) + 0.1)*ones(n,1);
t = 1;
for k=1:MAXITERS
    X = B+diag(x);  
    R = chol(X);  Rinv = inv(R);  Xinv = Rinv*Rinv';
    val  = t*sum(x) - 2*sum(log(diag(R))); 
    g = t - diag(Xinv);  
    H = Xinv.^2;
    dx = -H\g;  
    fprime = g'*dx;
    if ((-fprime/2) < NTTOL)
        Z = (1/t)*(Xinv - Xinv*diag(dx)*Xinv);
        gap = X(:)'*Z(:);    
        gaps = [gaps, gap];
        iters = [iters, k];  
        if (gap < tol);  return;  end; 
        t = min(t*mu, (n+1)/tol);  
    else
        s = 1;
        dX = Rinv'*diag(dx)*Rinv;
        Xnew = eye(n) + s*dX;  [Rnew,d] = chol(Xnew);
        while (d), 
            s = BETA*s; 
            Xnew = eye(n) + s*dX;
            [Rnew,d] = chol(Xnew);
        end;
        while (t*s*sum(dx) - 2*sum(log(diag(Rnew))) > ALPHA*s*fprime), 
            s = BETA*s; 
            Xnew = eye(n) + s*dX;   [Rnew,d] = chol(Xnew);
        end;
        x = x+s*dx;   
    end;
end;
