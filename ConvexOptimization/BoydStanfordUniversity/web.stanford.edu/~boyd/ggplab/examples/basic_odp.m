% Optimal doping profile optimization via GP
% (a figure is generated).
%
% This is an example taken directly from the paper:
%
%   Optimal Doping Profiles via Geometric Programming,
%   IEEE Transactions on Electron Devices, December, 2005,
%   by S. Joshi, S. Boyd, and R. Dutton.
%   (see pages 3-8 for problem details)
%
% Determines the optimal doping profile that minimizes base transit
% time in a (homojunction) bipolar junction transistor.
% This problem can be posed as a GP:
%
%   minimize   tau_B
%       s.t.   Nmin <= v <= Nmax
%              y_(i+1) + v_i^const1 <= y_i
%              w_(i+1) + v_i^const2 <= w_i, etc...
%
% where variables are v_i, y_i, and w_i.
%
% Almir Mutapcic and Siddharth Joshi 10/05
clear all;

% discretization size
M = 100;
% M = 1000; % takes a few minutes to process constraints

% problem constants
g1 = 0.42;
g2 = 0.69;
Nmax = 5*10^18;
Nmin = 5*10^16;
Nref = 10^17;
Dn0 = 20.72;
ni0= 1.4*(10^10);
WB = 10^(-5);
C =  WB^2/((M^2)*(Nref^g1)*Dn0);

% exponent powers
pwi = g2 -1;
pwj = 1+g1-g2;

% optimization variables
gpvar v(M) y(M) w(M) 

% problem constraints
constr = [ Nmin*ones(M,1) <= v; v <= Nmax*ones(M,1); ];

for i = 1:M-1
  if( mod(i,100) == 0 ), disp(i), end;
  constr(end+1) = y(i+1) + v(i)^pwj <= y(i);
  constr(end+1) = w(i+1) + y(i)*v(i)^pwi <= w(i);
end

constr(end+1) = y(M) == v(M)^pwj;
constr(end+1) = w(M) == y(M)*v(M)^pwi;

% objective function is the base transmit time
tau_B = C*w(1);

% solve the optimal doping profile problem
[opt_val sol status] = gpsolve(tau_B, constr);
assign(sol);

% plot the basic optimal doping profile
nbw = 0:1/M:1-1/M;
semilogy(nbw,v,'LineWidth',2);
axis([0 1 1e16 1e19]);
xlabel('base');
ylabel('doping');
text(0,Nmin,'Nmin ', 'HorizontalAlignment','right');
text(0,Nmax,'Nmax ', 'HorizontalAlignment','right');
