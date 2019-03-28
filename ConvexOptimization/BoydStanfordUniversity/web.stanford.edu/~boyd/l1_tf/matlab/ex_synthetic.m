%
% Comparison example of
%   l1 trend filtering and HP filtering on synthetic signal
%
%  first example given in "l1 Trend Filtering" Kim, Koh, Boyd and Gorinevsky

rand('state',72);
randn('state',72);

% data generation
n = 1000;       % total time
p = 0.01;       % change probability
sigma = 1;
noise = 10;

e = noise*randn(n,1);
v = sigma*(rand(n,1)-0.5);
for i = 2:n
    if (rand(1)>=p)
        v(i) = v(i-1);
    end
end
x = cumsum(v);
y = x + e;

% double difference matrix
I   = speye(n,n);
I2  = speye(n-2,n-2);
O2  = zeros(n-2,1);
D   = [I2 O2 O2]+[O2 -2*I2 O2]+[O2 O2 I2];

% regularization parameter
lambda = 1.25e+4;

%----------------------------------------------------------------------
%   l1 filter
%----------------------------------------------------------------------
[z1,status] = l1tf(y,lambda);

% uncomment line below to solve l1 trend filtering problem using CVX
% [z1,status] = l1tf_cvx(y,lambda);

err1 = norm(z1-y,2);

%----------------------------------------------------------------------
%   HP filter
%----------------------------------------------------------------------
L = 1e-1; U = 1e10;
for i = 1:100
    lambda = sqrt(L*U);
    if (lambda <= L || lambda >= U) break; end
    
    z2 = (speye(n)+lambda*D'*D)\y;
    
    err2 = norm(z2-y,2);
    if (err2 > err1)
        U = lambda;
    else
        L = lambda;
    end
end
disp(sprintf('error l1t = %e', err1));
disp(sprintf('error H-P = %e', err2));


%----------------------------------------------------------------------
%	Display result
%----------------------------------------------------------------------
tt   = 1:n;
xyzs = [x;y;z1;z2];
maxx = max(xyzs);
minx = min(xyzs);
SAVE_FIGURE = false; % set it true if you want save figures

if (SAVE_FIGURE)
    figure(1);
    plot(tt, x, '-'); ylim([minx maxx]);
    print -depsc synthetic_original
    figure(2);
    plot(tt,y); ylim([minx maxx]);
    print -depsc synthetic_noisy
    figure(3);
    plot(tt,z1,'r',tt,x,':k'); ylim([minx maxx]);
    print(3,'-depsc',sprintf('synthetic_l1t_%04d',floor(log10(lambda)*100)));
    figure(4);
    plot(tt,z2,'r',tt,x,':k'); ylim([minx maxx]);
    print(4,'-depsc',sprintf('synthetic_hp_%04d',floor(log10(lambda)*100)));
else
    figure(1);
    subplot(2,2,1); plot(tt,x,'-'); ylim([minx maxx]);
    subplot(2,2,2); plot(tt,y, 'b',tt,x,' k'); ylim([minx maxx]);
    subplot(2,2,3); plot(tt,z1,'r',tt,x,':k'); ylim([minx maxx]);
    subplot(2,2,4); plot(tt,z2,'r',tt,x,':k'); ylim([minx maxx]);
end
