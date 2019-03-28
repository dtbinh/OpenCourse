% Digital circuit sizing (a simple, hard-coded example).
% (a figure is generated if the tradeoff flag is turned on)
%
% This is an example taken directly from the paper:
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
% This code is specific to the digital circuit shown in figure 4
% (page 28) of GP tutorial paper. All the constraints and
% the worst-case delay expression are hard-coded for this
% particular circuit.
% 
% A more general code with more precise models for digital cicuit
% sizing is also available as part of the ggplab examples library.
% 
% Almir Mutapcic 02/01/2006
clear all; close all;
PLOT_TRADEOFF = 1; % to disable tradeoff plot, set PLOT_TRADEOFF = 0 

% number of cells
m = 7;

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
gpvar x(m)    % scale factors

% input capacitance is an affine function of sizes
cin = alpha + beta.*x;

% load capacitance of a gate is the sum of its fan-out c_in's
cload(1) = cin(4);
cload(2) = cin(4) + cin(5); 
cload(3) = cin(5) + cin(7); 
cload(4) = cin(6) + cin(7); 
cload(5) = cin(7); 
% output gates have their load capacitances
cload(6) = Cout6;
cload(7) = Cout7;

% gate delay is the product of its driving res. R = gamma./x and cload
d = (cload').*gamma./x;              

power = (f.*e)'*x;         % total power
area = a'*x;               % total area

% constraints
constr_x = ones(m,1) <= x; % all sizes greater than 1 (normalized)

% evaluate delay over all paths in the given circuit (there are 7 paths)
path_delays = [ ...
d(1) + d(4) + d(6); % delay of path 1
d(1) + d(4) + d(7); % delay of path 2, etc...
d(2) + d(4) + d(6);
d(2) + d(4) + d(7);
d(2) + d(5) + d(7);
d(3) + d(5) + d(6);
d(3) + d(7) ];

% objective is the worst-case delay
circuit_delay = max(path_delays); 

% collect all the constraints
constr = [power <= Pmax; area <= Amax; constr_x];

% formulate the GP problem and solve it
[obj_value, solution, status] = gpsolve(circuit_delay, constr);
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
N = 25;
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
    [obj_value, solution, status] = gpsolve(circuit_delay, constr);
    min_delay(k,n) = obj_value;
  end
end

% enable solver reporting again
global QUIET; QUIET = 0;

% plot the tradeoff curve
plot(Pmax,min_delay(1,:), Pmax,min_delay(2,:), Pmax,min_delay(3,:));
xlabel('Pmax'); ylabel('Dmin');

end
% end tradeoff curve code
