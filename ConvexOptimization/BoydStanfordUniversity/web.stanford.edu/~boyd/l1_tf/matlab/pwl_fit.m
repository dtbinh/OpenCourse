function x = pwl_fit(y,idx)
% x = pwl_fit(y,idx)
%
% finds a piecewise fit for y given change point index
%
% INPUT:
%   y       : n-vector; original signal
%   idx     : k-vector; index of change point
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

x = P*((P'*P)\(P'*y));
