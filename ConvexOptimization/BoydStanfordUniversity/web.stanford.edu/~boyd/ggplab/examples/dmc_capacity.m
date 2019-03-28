% Discrete memoryless channel (DMC) capacity computation.
%
% For more details see:
%
%    Geometric Programming Duals of Channel Capacity and Rate Distortion
%    by M. Chiang and S. Boyd,
%    IEEE Transactions on Information Theory, 2004.
%
% Computes exact capacity of a discrete memoryless channel (DMC)
% via the dual problem, which is a GP.  GGPLAB does not return
% dual variables, so we can't retrieve the optimal input distribution.
% The Lagrange dual of the channel capacity problem is GP:
%
%   minimize   sum(z)
%       s.t.   prod(z_j^{P_ij}) >= exp(-entr(P^(i))),  for i = 1,...,N
%
% where variable is z, and P^(i) is the ith row of the channel matrix P.
%
% Almir Mutapcic 01/18/2005

% transition probability matrix for a discrete memoryless channel
P = [0.3 0.2 0.5; 0.5 0.3 0.2; 0.2 0.5 0.3];

% number of input (N) and output (M) symbols
[N,M] = size(P);

% GP variables
gpvar z(M)

% objective is the channel capacity
obj = sum(z);

% constraints
constr = gpconstraint;
for k = 1:N
  constr(k) = exp( -entropy(P(k,:)) ) <= prod( z.^(P(k,:)') );
end

% find the channel capacity
[capacity sol status] = gpsolve(obj,constr);

disp(' ')
disp(['Channel capacity (via dual) is ' num2str(capacity) '.'])
