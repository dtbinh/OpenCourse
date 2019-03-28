function x = pwl_fit_l1(y,idx,x0,lambda)
% x = pwl_fit_l1(y,idx,x0,lambda)
%
% finds a piecewise fit for y given change point index
%
% INPUT:
%   y       : n-vector; original signal
%   idx     : k-vector; index of change point
%   x0      : n-vector; approximate solution of the l1 trend estimation
%   lambda  : scalar; regularization parameter
%
% OUTPUT:
%   x       : n-vector; piecewise linear fit of y
%

n   = length(y);      % length of signal x
I2  = speye(n-2,n-2);
O2  = zeros(n-2,1);
D   = [I2 O2 O2]+[O2 -2*I2 O2]+[O2 O2 I2];

% add both end points if not included
if (idx(1) ~= 1)  , idx = [1; idx]; end
if (idx(end) ~= n), idx = [idx; n]; end

r = length(idx);
P = zeros(n,r);
for i = 1:r-1
    prev = idx(i);
    next = idx(i+1);
    k = [prev:next-1]';
    P(k,i:i+1)   = [(next-k)/(next-prev),(k-prev)/(next-prev)];
end
P(end,end) = 1;
P = sparse(P);

dist= diff(idx);
d1  = 1./dist(1:end-1);
d2  = 1./dist(2:end);
v   = x0(idx);

G   = spdiags([d1 -(d1+d2) d2],0:2,r-2,r);

vv  = ((P'*P)\(P'*y-lambda*(G'*sign(G*v))));

if ( sign(G*v)~=sign(G*vv) )
    disp('possible sign change in v');
    disp('try l1t_primal_trim');
end

x = P*vv;
