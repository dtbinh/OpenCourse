% Optimal two-pole lowpass filter implementation.
% (a figure with the tradeoff between the power and noise is generated)
% FIX: make sure that this numerical example is fine.
%
% This is an example taken directly from the paper:
%
%   Geometric Programming and its Applications to EDA Problems
%   (DATE Tutorial 2005) by Boyd, Kim, and Mohan.
%   (see pages 87-92)
%
% Implements an optimal two-pole Butterworth filter with frequency wc,
% which minimizes total power, subject to area, output noise limits, etc.
% The problem can be posed as a generalized GP:
%
%   minimize   P(u1) + P(u2) 
%       s.t.   t1 = sqrt(2)/wc,   t2 = (1/sqrt(2))/wc
%              Aamp(u1) + Aamp(u2) + Acap(v1) + Acap(v2) <= Amax
%              M = (wc/4*sqrt(2))*sqrt(N1^2 + 2*N2^2) <= Mmax
%
% where variables are u1, u2, v1, and v2.
%
% Almir Mutapcic 02/02/06
close all; clear all;

% set the quiet flag (no solver reporting)
global QUIET; QUIET = 1;

% problem constants
wc = 10^8;
wc = 100;
Amax = 4*10^-6;
Amax = 4*10^-1;

% GP variables
gpvar W1 L1 W2 L2 Acap1 Acap2

% amplifier model
Aamp1 = W1*L1;
Aamp2 = W2*L2;

P1 = (2.5*10^-4)*(W1/L1);
P2 = (2.5*10^-4)*(W2/L2);

g1 = (4*10^-5)*W1/L1;
g2 = (4*10^-5)*W2/L2;

N1 = sqrt(7.5*10^-16)*sqrt(L1/W1);
N2 = sqrt(7.5*10^-16)*sqrt(L2/W2);

C1 = (10^-4)*Acap1;
C2 = (10^-4)*Acap2;

t1 = C1/g1;
t2 = C2/g2;

M = wc/(4*sqrt(2))*sqrt(N1^2 + 2*N2^2);

% variable parameter (max noise)
Mmax = [10:10:100]*10^-6;
Mmax = [10:10:100]*10^-4;

% objective is to minimize total power of the implementation
% obj = P(u1) + P(u2) (private variables, let u1 = L and u2 = W)
obj = P1 + P2;

% constraint set 
constr = [ ...
  t1 == sqrt(2)/wc;
  t2 == (1/sqrt(2))/wc;
  Aamp1 + Aamp2 + Acap1 + Acap2 <= Amax;
];

power = zeros(1,10);
status = {};
% setup GP problems with varying parameter Mmax and obtain tradeoff curve
for m = 1:length(Mmax)
  % add constraints that have varying parameters
  constr(4) = M <= Mmax(m);

  % solve the optimal filter implementation problem
  prob = gpproblem(obj, constr);
  res = solve(prob);
  power(m) = res.obj_value;
  status{end+1} = res.status;
end

disp('Displaying power vector: ')
power

% enable solver reports
global QUIET; QUIET = 0;

% tradeoff plot
semilogy(Mmax*10^3, power*10^3);
xlabel('max noise Mmax (milli-Volts)');
ylabel('power P (milli-Watts)');
