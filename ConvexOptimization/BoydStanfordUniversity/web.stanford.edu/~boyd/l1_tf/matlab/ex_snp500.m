%
% Comparison example of 
%   l1 trend filtering and HP filtering on S&P500 data
%
%  second example given in "l1 Trend Filtering" Kim, Koh, Boyd and Gorinevsky


% load data
data    = csvread('snp500.dat');
y       = log(data(end:-1:1,end));  % reverse order, last column, log-scale
dates   = data(end:-1:1,1:3);       % reverse order, 1~3 column
n       = length(y);

% regularization parameter
lambda  = 50;

I       = speye(n,n);
I2      = speye(n-2,n-2);
O2      = zeros(n-2,1);
D       = [I2 O2 O2]+[O2 -2*I2 O2]+[O2 O2 I2];

%----------------------------------------------------------------------
%   l1 trend filter
%----------------------------------------------------------------------
[z1] = l1tf(y,lambda);

% uncomment line below to solve l1 trend filtering problem using CVX
% [z1] = l1tf_cvx(y,lambda);

err1 = norm(z1-y,2);

%----------------------------------------------------------------------
%   HP filter
%----------------------------------------------------------------------
% solve HP filtering problem with lambda that gives the same fitting
% error in l2-norm with l1 trend filtering, using bisection.

L = 1e-1; U = 1e10;
for i = 1:100
    lambda2 = sqrt(L*U);
    if (lambda2 <= L || lambda2 >= U) break; end
    
    z2 = (speye(n)+lambda2*D'*D)\y;
    
    err2 = norm(z2-y,2);
    
    if (err2 > err1)
        U = lambda2;
    else
        L = lambda2;
    end
end

%----------------------------------------------------------------------
%   Display result
%----------------------------------------------------------------------
disp(sprintf('error l1t = %e', err1));
disp(sprintf('error H-P = %e', err2));

plot_signals(1:n,y,z1,z2,'snp500',lambda,1,dates);
