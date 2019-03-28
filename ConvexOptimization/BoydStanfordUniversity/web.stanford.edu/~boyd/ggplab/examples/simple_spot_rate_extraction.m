% Simplified spot rate curve extraction (coupon stripping).
%
% This is a simplified example from the paper:
%
%   Analyzing the Indonesian Government Bond Market
%   with a New Coupon Stripping Model and Principal Components
%   by K. O. Kortanek and H. Yunianto
%   (see page 15 for the GP problem formulation)
%
% Performs the spot rate curve extraction (coupon stripping) from
% financial data that only include bonds with no coupon payments.
% This results in a GP that can be found on the page 15 of the paper.
%
% Almir Mutapcic 01/15/06

%********************************************************************
% problem data
%********************************************************************
% problem sizes
N = 2;   % number of bonds (length of PVhat)
L = 1;   % number of uniform perturbation intervals

% observed prices based on hypothetical bond data (this is made up data)
PVhat = [98.85; 96.25];
% maturity for each bond (in years)
T = [0.008; 0.015];

% problem constants (all made up)
RH    = 1; RL    = .5;
RHinf = 1; RLinf = .5;
WH    = 1; WL    = .5;

global beta gamma1 gamma2;
beta   = -0.1005;
gamma1 = 0.999*ones(N,1);
gamma2 = 0.999*ones(N,1);

%********************************************************************
% problem variables, objective function, and constraints
%********************************************************************
% GP variables
gpvar x(2*L+2+1)     % note: w_bar = exp(x_2*L+2+1)

% objective function
obj = inv(x(L+2+1)) + x(L+2+1);

% model prices
% denote P_k(T_k) as PV
PV = posynomial;
for k = 1:N
  PV(k,1) = 100*x(1)^ar(T(k))*x(2)^aalpha(T(k));
  PV(k,1) = PV(k)*x(3)^b(T(k)); % here b is the same as aalpha
end

% hard constraint on r
c1 = [ exp(-RH)*x(1) <= 1; exp(RL)*x(1)^(-1) <= 1];

% hard constraint on alpha/-beta
c2 = [ exp(RHinf*beta)*x(2) <= 1; exp(-RLinf*beta)*x(2)^(-1) <= 1];

% "driving" constraints
c_drive = [ exp(-WH)*x(3)*x(4)^-1*x(5)^-1 <= 1;
            exp(WL)*x(3)^-1*x(4)^-1*x(5)  <= 1 ];

% constraint set 
constr = [ c_drive; c1; c2; PV <= PVhat ];

%********************************************************************
% solving the GP problem
%********************************************************************
% extract the spot rate curve from the observables
[oval sol status] = gpsolve(obj, constr);

% evaluate and display computed prices
disp(' ')
fprintf(1,'Estimate of PV(1) is %3.4f versus observed PV_hat(1) = %3.4f.\n',...
        eval(PV(1), sol), PVhat(1));
fprintf(1,'Estimate of PV(2) is %3.4f versus observed PV_hat(2) = %3.4f.\n',...
        eval(PV(2), sol), PVhat(2));
