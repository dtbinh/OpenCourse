% Digital circuit sizing using vectorize features
% (a figure is generated if the tradeoff flag is turned on)
%
% This is an simple, ahard-coded example taken directly from:
%
%   A Tutorial on Geometric Programming (see pages 25-29)
%   by Boyd, Kim, Vandenberghe, and Hassibi.
%
% Solves the problem of choosing gate scale factors x_i to give
% minimum ckt delay, subject to limits on the total area and power.
%
%   minimize   D
%       s.t.   P <= Pmax, A <= Amax
%              x >= 1
%
% where variables are scale factors x.
%
% This code uses matrices in order to evaluate signal paths
% through the circuit (thus, it uses vectorize Matlab features).
% It is specific to the digital circuit shown in figure 4 (page 28)
% of GP tutorial paper.
%
% Almir Mutapcic 02/01/2006
clear all; close all;
PLOT_TRADEOFF = 1; % to disable set PLOT_TRADEOFF = 0;

% digital circuit shown in figure 4 (page 28) of GP tutorial paper
m = 7;  % number of cells
n = 8;  % number of edges
A = sparse(m,n);

% A is standard cell-edge incidence matrix of the circuit
% A_ij = 1 if edge j comes out of cell i, -1 if it comes in, 0 otherwise
  A(1,1) =     1;
  A(2,2) =     1;
  A(2,3) =     1;
  A(3,4) =     1;
  A(3,8) =     1;
  A(4,1) =    -1;
  A(4,2) =    -1;
  A(4,5) =     1;
  A(4,6) =     1;
  A(5,3) =    -1;
  A(5,4) =    -1;
  A(5,7) =     1;
  A(6,5) =    -1;
  A(7,6) =    -1;
  A(7,7) =    -1;
  A(7,8) =    -1;

% decompose A into edge outgoing and edge-incoming part 
Aout = double(A > 0);
Ain = double(A < 0);

% problem constants
f = [1 0.8 1 0.7 0.7 0.5 0.5]';
e = [1 2 1 1.5 1.5 1 2]';
Cout6 = 10;
Cout7 = 10;

a     = ones(m,1);
alpha = ones(m,1);
beta  = ones(m,1);
gamma = ones(m,1);

% maximum area and power specification 
Amax = 25; Pmax = 50;

% optimization variables
gpvar x(m)                 % sizes
gpvar t(m)                 % arrival times

% input capacitance is an affine function of sizes
cin = alpha + beta.*x;

% load capacitance is the input capacitance times the fan-out matrix
% given by Fout = Aout*Ain' 
cload = (Aout*Ain')*cin;
cload(6) = Cout6;          % load capacitance of the output gate 6
cload(7) = Cout7;          % load capacitance of othe utput gate 7

% delay is the product of its driving resistance R = gamma./x and cload
d = cload.*gamma./x;              

power = (f.*e)'*x;         % total power
area = a'*x;               % total area

% constraints
constr_x = ones(m,1) <= x; % all sizes greater than 1 (normalized)

% create timing constraints
% these constraints enforce t_j + d_j <= t_i over all gates j that drive gate i
constr_timing = Aout'*t + Ain'*d <= Ain'*t;
% for gates with inputs not connected to other gates we enforce d_i <= t_i
input_timing  = d(1:3) <= t(1:3);

% objective is the upper bound on the overall delay 
% and that is the max of arrival times for output gates 6 and 7
D = max(t(6),t(7));

% collect all the constraints
constr = [power <= Pmax; area <= Amax; constr_timing; input_timing; constr_x];

% formulate the GP problem and solve it
[obj_value, solution, status] = gpsolve(D, constr);
assign(solution);

fprintf(1,'\nOptimal circuit delay for Pmax = %d and Amax = %d is %3.4f.\n', ...
        Pmax, Amax, obj_value)
disp('Optimal scale factors are: ')
x

%********************************************************************
% tradeoff curve code
%********************************************************************
if( PLOT_TRADEOFF )

% varying parameters for an optimal trade-off curve
N = 10;
Pmax = linspace(10,100,N);
Amax = [25 50 100];
min_delay = zeros(length(Amax),N);

% set the quiet flag (no solver reporting)
global QUIET; QUIET = 1;

for k = 1:length(Amax)
  for n = 1:N
    % add constraints that have varying parameters
    constr(1) = power <= Pmax(n);
    constr(2) = area  <= Amax(k);

    % solve the GP problem and compute the optimal volume
    [obj_value, solution, status] = gpsolve(D, constr);
    min_delay(k,n) = obj_value;
  end
end

% enable solver reporting again
global QUIET; QUIET = 0;

plot(Pmax,min_delay(1,:), Pmax,min_delay(2,:), Pmax,min_delay(3,:));
xlabel('Pmax'); ylabel('Dmin');

end
% end tradeoff curve code
