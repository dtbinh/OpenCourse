function [L, d, flag, iter] = toeplitzFact(r,pd)
% [L,d,flag, iter] = TOEPLITZFACT(r,pd)
% computes the factorization of a nonsingular symmetric Toeplitz matrix
% L * toeplitz(r) * L' = diag(d)
% where L is unit lower triangular and D is diagonal. 
% The flag indicates whether toeplitz(r) is positive definite (flag = 1)
% or not (flag = 0).
%
% Procedure stops early and outputs L = 0 if toeplitz(r) is singular,
% or is not positive definite when required. The number of iterations run
% is iter.
%
% Author: Hsiao-Han (Suzi) Chao (suzihchao@ucla.edu)
% Date modified: 7/18/2018

n = length(r);
if nargin == 1
    pd = false; % does not require positive definiteness
end
d = zeros(n,1);
L = eye(n);
flag = 1; % is positive definite
if r(1) <= 0
    flag = 0; 
    if pd || (r(1) == 0), L = 0; iter = 1; return; end
end
d(1) = r(1);
for k = 2:n
    delta = L(k-1,1:k-1)*r(2:k);
    d(k) = d(k-1)-delta^2/d(k-1);
    if d(k) <= 0
        flag = 0; 
        if pd || (d(k) == 0), L = 0; iter = k; return; end
    end
    L(k,k-1:-1:1) = [L(k-1,k-2:-1:1) 0] - (delta/d(k-1))*L(k-1,1:k-1);
end
iter = n;
end