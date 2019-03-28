% Optimization of inductor circuits via geometric programming (GP).
%
% TODO: This example is not complete, it still does not reproduce
%       plots shown in the reference paper.
% FIX:  The problem is primal infeasible as reported by Mosek.
%       We are missing some constants necessary for the solution,
%       which were not fully specified in the paper.
%
% This example is taken directly from the paper:
%
%   Optimization of inductor circuits via geometric programming
%   Proceedings of the 36th {ACM/IEEE} Conference on Design Automation
%   New Orleans, June 1999 (pages 994-998)
%   by M. Hershenson, S. Mohan, S. Boyd, and T. Lee 
%   (see pages 995-996 for problem details)
%
% Designs a sprial inductor that maximizes the quality factor for a given
% inductance value and for given constraints on the maximum area, minimum
% spacing and turn width, etc. The problem can be posed as a GP:
%
%   maximize   Q_L,min
%       s.t.   Q_L >= Q_L,min
%              L == Lreq, etc...
%
% where variables are number of turns n, the turn width w, the turn spacing s,
% the outer diameter d_out, the average diameter d_avg, and min quality factor.
%
% Almir Mutapcic 12/05

% problem constants
A_min = 400; % 400 micro m^2
w_min= 1.9;
s_min = 1.9;
freq = 2.5; % 2.5 GHz
w_sr_min = 7;

% optimization variables
gpvar n w s d_out d_avg QL_min

% inductor area
A = d_out^2;

% inductor length
l = 4*n*d_avg;

% series resistance Rs
sigma = 3*10^5;
t_M5 = 0.9;
mu = 4*pi*10^-7;
%%% delta = sqrt(2/(freq*mu*sigma));        % skin depth
%%% k1 = 1/(sigma*delta*(1-exp(-t_M5/delta)));
k1 = 0.001;
Rs = k1*l/w;

% series inductance
beta = 1.66*10^-3;
alpha1 = -1.33;
alpha2 = -0.125;
alpha3 = 2.5;
alpha4 = 1.83;
alpha5 = -0.022;
Ls = beta*d_out^(alpha1)*w^(alpha2)*d_avg^(alpha3)*n^(alpha4)*s^(alpha5);

% minimum self-resonance constraint
% (will consider the unconstrained case)

% quality factor constraint
QL = freq*Ls/Rs;

% problem constraints
constr = [ w_min <= w; s_min <= s; A_min <= A; QL_min <= QL ];

% maximum constraints (adding on my own)
% constr = [ constr; w <= 2*w_min; s <= 2*s_min; A <= 2*A_min ];

% create the changing constraint
L_req = 0;
constr(end+1) = Ls == L_req;

% sweep the inductance L_req
Qarray = [];
% for L_req = 2:0.5:20
for L_req = 2:0.5:3
  constr(end) = Ls == L_req;
  show(constr)
  % maximize the minimim quality factor QL_min 
  [opt_val sol status] = gpsolve(QL_min, constr, 'max');
  Qarray = [Qarray opt_val];
end

% FIX: will also need to round to quarter sizes

% plot the optimal trade-off curve
% L_req = 2:0.5:20;
% plot(L_req,Qarray);
