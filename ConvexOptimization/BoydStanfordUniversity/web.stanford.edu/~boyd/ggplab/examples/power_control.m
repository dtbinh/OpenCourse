% Simple power control in communication systems via GP.
%
% This is an example taken from the GP tutorial paper:
%
%   A Tutorial on Geometric Programming (see pages 16-17)
%   by Boyd, Kim, Vandenberghe, and Hassibi.
%
% Solves the power control problem in communication systems, where
% we want to minimize the total transmitter power for n transmitters,
% subject to minimum SINR level, and lower and upper bounds on powers.
% This results in a GP: 
%
%   minimize   sum(P)
%       s.t.   Pmin <= P <= Pmax
%              SINR >= SINR_min
%
% where variables are transmitter powers P.
% Numerical data for the specific examples was made up.
%
% Almir Mutapcic 01/15/06

% problem constants
n = 5;                 % number of transmitters and receivers
sigma = 0.5*ones(n,1); % noise power at the receiver i
Pmin = 0.1*ones(n,1);  % minimum power at the transmitter i
Pmax = 5*ones(n,1);    % maximum power at the transmitter i
SINR_min = 2;          % threshold SINR for each receiver

% path gain matrix
G = [1.0  0.1  0.2  0.1  0.0
     0.1  1.0  0.1  0.1  0.0
     0.2  0.1  2.0  0.2  0.2
     0.1  0.1  0.2  1.0  0.1
     0.0  0.0  0.2  0.1  1.0];

% variables are power levels
gpvar P(n) 

% objective function is the total transmitter power
Ptotal = sum(P);

% formulate the inverse SINR at each receiver using vectorize features
Gdiag = diag(G);          % the main diagonal of G matrix
Gtilde = G - diag(Gdiag); % G matrix without the main diagonal
% inverse SINR
inverseSINR = (sigma + Gtilde*P)./(Gdiag.*P);

% constraints are power limits and minimum SINR level 
constr = [ Pmin <= P; P <= Pmax;
           inverseSINR <= (1/SINR_min)*ones(n,1); ];

% solve the power control problem
[min_Ptotal solution status] = gpsolve(Ptotal, constr);
assign(solution);

fprintf(1,'\nThe minimum total transmitter power is %3.2f.\n',min_Ptotal);
disp('Optimal power levels are: '), P
