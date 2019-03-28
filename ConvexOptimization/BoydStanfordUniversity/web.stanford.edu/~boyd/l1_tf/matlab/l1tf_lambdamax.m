function lambdamax = l1tf_lambdamax(y)
% lambdamax = l1tf_lambdamax(y)
%
% returns an upperbound of lambda:
%   With a regularization parameter value over lambda_max, l1tf returns the
%   best affine fit for y.
% 
% INPUT:
%   y           : n-vector; original signal
%
% OUTPUT:
%
%   lambdamax   : scalar; maximum value of lambda in useful range
%

n   = length(y);
I   = speye(n-2,n-2);
O   = zeros(n-2,1);
D   = [I O O]+[O -2*I O]+[O O I];

lambdamax = norm((D*D')\(D*y),inf);

disp(sprintf('  lambda_max : %e', lambdamax)) ;
