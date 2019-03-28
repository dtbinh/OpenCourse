% Logistic regression modeling via geometric programming (GP).
% (a figure is generated)
%
% This examples solves a logistic regression example presented
% in the book "Convex Optimization" by Boyd and Vandenberghe
% (see pages 354-355). More info can be found in the attached report:
%
%   Logistic regression via Geometric Programming
%   by Seung Jean Kim and Almir Mutapcic
%   (Will be available soon.)
%
% Solves the logistic regression problem re-formulated as a GP.
% The original log regression problem is:
%
%   minimize   sum_i(theta'*x_i) + sum_i( log(1 + exp(-theta'*x_i)) )
%
% where x are explanatory variables and theta are model parameters.
% The equivalent GP is obtained by the following change of variables:
% z_i = exp(theta_i). The log regression problem is then a GP:
%
%   minimize   prod( prod(z_j^x_j) ) * (prod( 1 + prod(z_j^(-x_j)) ))
%
% with variables z and data x (explanatory variables).
%
% Almir Mutapcic, 11/05

% load problem data from the Convex Optimization book
load_log_reg_data;

% order the observation data
ind_false = find( y == 0 );
ind_true  = find( y == 1 );

% X is the sorted design matrix
% first have true than false observations followed by the bias term
X = [u(ind_true); u(ind_false)];
X = [X ones(size(u,1),1)];
[m,n] = size(X);
q = length(ind_true);

% optimization variables
gpvar z(n) t(q) s(m)

% objective function
obj = prod(t)*prod(s);

constr = gpconstraint;
% constraints
for k = 1:q
  constr(k) = prod( z.^(X(k,:)') ) <= t(k);
end

for k = 1:m
  constr(end+1) = 1 + prod( z.^(-X(k,:)') ) <= s(k);
end

% solve the GP problem
[obj_value, solution, status] = gpsolve(obj, constr)
assign(solution)

% retrieve the optimal values and plot the result
theta = log(z);
aml = -theta(1);
bml = -theta(2);

us = linspace(-1,11,1000)';
ps = exp(aml*us + bml)./(1+exp(aml*us+bml));

plot(us,ps,'-', u(ind_true),y(ind_true),'o', ...
                u(ind_false),y(ind_false),'o');
axis([-1, 11,-0.1,1.1]);
