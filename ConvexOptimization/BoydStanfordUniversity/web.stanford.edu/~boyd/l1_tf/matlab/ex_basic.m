%
%   Simple usage example of l1tf
%
rand('state',0); randn('state',0); % make data reproducible

%----------------------------------------------------------------------
% 	data generation
%----------------------------------------------------------------------
x = interp1([1 100 150 260 300 350 400],[3 10 -2 -5 9 2 6],1:400)';
n = length(x);
m = n-2;
y = x+3*randn(n,1);

% set maximum regularization parameter lambda_max
lambda_max = l1tf_lambdamax(y);

%----------------------------------------------------------------------
% 	l1 trend filtering
%----------------------------------------------------------------------
[z1,status] = l1tf(y, 0.001*lambda_max);
[z2,status] = l1tf(y, 0.005*lambda_max);
[z3,status] = l1tf(y, 0.025*lambda_max);

% uncomment line below to solve l1 trend filtering problem using CVX
% [z1,status] = l1tf_cvx(y, 0.001*lambda_max);
% [z2,status] = l1tf_cvx(y, 0.005*lambda_max);
% [z3,status] = l1tf_cvx(y, 0.025*lambda_max);

%----------------------------------------------------------------------
% 	plot results
%----------------------------------------------------------------------
xyzs = [x y z1 z2 z3];
maxx = max(max(xyzs));
minx = min(min(xyzs));

figure(1);
subplot(2,2,1); plot(1:n,x,'b',1:n,y, 'r'); ylim([minx maxx]); title('original');
subplot(2,2,2); plot(1:n,x,'b',1:n,z1,'r'); ylim([minx maxx]); title('lambda = 0.001');
subplot(2,2,3); plot(1:n,x,'b',1:n,z2,'r'); ylim([minx maxx]); title('lambda = 0.005');
subplot(2,2,4); plot(1:n,x,'b',1:n,z3,'r'); ylim([minx maxx]); title('lambda = 0.025');

