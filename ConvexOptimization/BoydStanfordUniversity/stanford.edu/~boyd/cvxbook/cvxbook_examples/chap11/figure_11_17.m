% Generates figures 11.17 and 11.18 of Boyd & Vandenberghe, Convex 
% Optimization
%
% Barrier method for a small SDP.

function figure_11_17()

% Generate a pair of strictly primal and dual feasible SDPs
%
% (primal) minimize    c'*x
%          subject to  x1*A1 + ... + xn*An <= B
%
% (dual)   maximize    -tr(B*Z)
%          subject to  tr(Ai*Z) + ci = 0, i=1,...,m
%                      Z >= 0

randn('state',0);
rand('state',0);
m = 100;
n = 50;

% Random symmetric Ai's 
A = randn(m^2,n);
for i=1:n
    mtrx = reshape(A(:,i),m,m);
    mtrx = (mtrx + mtrx')/2;
    A(:,i) = mtrx(:);
end;

% Random positive definite B
B = randn(m,m);  B = B*B';  [V,D] = eig(B);  B = V*diag(rand(m,1))*V';

% ci = -tr(Ai*Z) for a random positive definite Z
Z = randn(m,m);  Z = Z*Z';  [V,D] = eig(Z);  Z = V*diag(rand(m,1))*V';
c = -A'*Z(:);

% Compute optimal vaue and shift B so that optimal value is 1.
x0 = zeros(n,1);
[x, Z, inniters, gaps] = sdp(A, B, c, x0, 10, 1e-6);
s = (1-c'*x)/norm(c)^2;
B = B + s*reshape(A*c,m,m);;
x0 = s*c; 


% Figure 11.17.  Gap versus Newton iterations for three values of mu.

figure(1)
muvals = [2, 50, 150];
[x, Z, iters, gaps] = sdp(A, B, c, x0, muvals(1),1e-6);
l = length(gaps);  iters1 = [];  gaps1 = [];
for i=1:l-1
    iters1 = [iters1, iters(i)-1, iters(i+1)-1];  
    gaps1 = [gaps1, gaps(i), gaps(i)]; 
end;
iters1 = [iters1, iters(l)-1] - iters1(1);  
gaps1 = [gaps1, gaps(l)]; 

[x, Z, iters, gaps] = sdp(A, B, c, x0, muvals(2), 1e-6);
l = length(gaps);  iters2 = [];  gaps2 = [];
for i=1:l-1
    iters2 = [iters2, iters(i)-1, iters(i+1)-1];  
    gaps2 = [gaps2, gaps(i), gaps(i)]; 
end;
iters2 = [iters2, iters(l)-1] - iters2(1);  
gaps2 = [gaps2, gaps(l)]; 

[x, Z, iters, gaps] = sdp(A, B, c, x0, muvals(3), 1e-6);
l = length(gaps);  iters3 = [];  gaps3 = [];
for i=1:l-1
    iters3 = [iters3, iters(i)-1, iters(i+1)-1];  
    gaps3 = [gaps3, gaps(i), gaps(i)]; 
end;
iters3 = [iters3, iters(l)-1] - iters3(1);  
gaps3 = [gaps3, gaps(l)]; 

semilogy(iters1,gaps1, iters2,gaps2, '-', iters3,gaps3, '--');
axis([0 100 0.5e-7 1e3]);
xlabel('Newton iterations'); ylabel('duality gap');
text(iters1(length(iters1)), gaps1(length(iters1)),'mu=2');
text(iters2(length(iters2)), gaps2(length(iters2)),'mu=50');
text(iters3(length(iters3)), gaps3(length(iters3)),'mu=150');


% Figure 11.19.  Number of Newton iterations versus mu.

figure(2)
muvals = [1.5 2 4 6 8 10:10:120];
noiters=zeros(1,length(muvals));
for i=1:length(muvals)
    [x, Z, iters, gaps] = sdp(A, B, c, x0, muvals(i), 1e-3);
    noiters(i) = iters(length(iters))-iters(1); % skip initial centering
    disp(['mu = ', num2str(muvals(i)), ': ', num2str(noiters(i)), ...
        ' Newton iterations.'])
end;
plot(muvals,noiters,'o', muvals,noiters,'-');
axis([ 0 120 0 140]);
xlabel('mu'); ylabel('Newton iterations');


function [x, Z, inniters, gaps] = sdp(A, B, c, x0, mu, tol)

% [x, Z, inniters, gaps] = sdp(A, B, c, x0, mu, tol)
%
% Solves the SDP 
% 
%     minimize  c'*x 
%     subjec to x1*A1 + ... + xn*An <= B
%
% x0:  strictly feasible starting point, not necesarily on central path
% B:   symmetric nxn matrix
% A:   [A1(:), A2(:) ... An(:)]
% inniters:  no of Newton iters per outer iteration
% gaps:   duality gap at the end of each outer iterations
 

n = length(c);
m = size(B,1);

MAXITERS = 500;   
ALPHA = 0.01;
BETA = 0.5;
NTTOL = 1e-5;     % stop inner iteration if lambda^2/2 < NTTOL

t = 1;
x = x0;
gaps = [];
inniters = [];
Asc = zeros(m^2,1);
for k=1:MAXITERS
    D = B - reshape(A*x,m,m);   
    R = chol(D);  Rinv = inv(R);  %D = R'*R 
    val = t*c'*x - 2*sum(log(diag(R)));
    for i=1:n
        Asc(:,i) = reshape(Rinv'*reshape(A(:,i),m,m)*Rinv, m^2,1);
    end;
    g = t*c + Asc'*reshape(eye(m),m^2,1);
    H = Asc'*Asc;
    dx = -H\g;  
    fprime = g'*dx;
    if ((-fprime/2) < NTTOL)
        Z = Rinv*Rinv'/t;
        gap = Z(:)'*D(:);    
        gaps = [gaps,gap];
        inniters = [inniters, k];  
        if (gap < tol), return; end; 
        t = min(t*mu, (m+1)/tol);  
    else
        s = 1;
        [Rnew, nonpd] = chol(B- reshape(A*(x+s*dx),m,m));
        while (nonpd), 
            s = BETA*s;
            [Rnew,nonpd] = chol(B- reshape(A*(x+s*dx),m,m));
        end;
        while (t*c'*(x+s*dx) - 2*sum(log(diag(Rnew))) > ...
            val + ALPHA*s*fprime)
            s = BETA*s; 
            Rnew = chol(B- reshape(A*(x+s*dx),m,m));
        end;
        x = x+s*dx;
    end;
end;
