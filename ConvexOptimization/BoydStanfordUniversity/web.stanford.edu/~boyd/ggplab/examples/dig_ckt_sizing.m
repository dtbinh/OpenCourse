% Digital circuit sizing example.
% (a figure is generated if the tradeoff flag is turned on)
%
% This is an example taken directly from the paper:
%
%   Digital Circuit Optimization via Geometric Programming
%   by Boyd, Kim, Patil, and Horowitz.
%   Operations Research 53(6): 899-932, 2005.
%
% Solves the problem of choosing gate scale factors x_i to give
% minimum ckt delay, subject to limits on the total area and power.
% Uses max gate arrival time T formulation that avoids evaluation
% of the delay over all paths in the circuit.
%
%   minimize   T_bar
%       s.t.   T_j <= T_bar      for j an output gate
%              T_j + d_i <= T_i  for j in FI(i)
%              P <= Pmax, A <= Amax
%              x >= 1
%
% where variables are x and T.
% 
% We use the circuit topology presented in figure 1 (page 902),
% where we take gates 1, 3 and 6 to be inverters (INV), 
% gates 2 and 7 to be three input NANDs (NAND3),
% and gates 4 and 5 to be two input NORs (NOR2).
%
% Almir Mutapcic 01/26/06
clear all; close all;
PLOT_TRADEOFF = 1; % to disable plotting set PLOT_TRADEOFF = 0;

%********************************************************************
% user specified data (specify problem constant and ckt topology)
%********************************************************************
m = 7;        % number of gates
Vdd = 5;      % supply voltage
Amax = 250;   % maximum area spec
Pmax = 15;    % maximum power spec

% gate specs
INV   = struct('Cin',3, 'Cint',3, 'Rdrv',0.48, 'A',3,  'Ileak',0.006);
NAND3 = struct('Cin',4, 'Cint',6, 'Rdrv',0.48, 'A',8,  'Ileak',0.007);
NOR2  = struct('Cin',5, 'Cint',6, 'Rdrv',0.48, 'A',10, 'Ileak',0.009);

gates([1 3 6]) = INV;
gates([2 7])   = NAND3;
gates([4 5])   = NOR2;

% primary inputs and primary outputs labels (start with m+1)
primary_inputs = [8 9 10];
primary_outputs = [11 12];
M = m + length( primary_inputs ) + length( primary_outputs );

% fan-in cell array
FI = cell(M,1);
FI{1} = [8];
FI{2} = [8 9 10];
FI{3} = [10];
FI{4} = [1 2];
FI{5} = [2 3];
FI{6} = [4];
FI{7} = [3 4 5];
FI{8} = [];
FI{9} = [];
FI{10} = [];
FI{11} = [6];
FI{12} = [7];

% primary output has Cin capacitance (but has no Cload)
Cin_po = sparse(M,1);
Cin_po(primary_outputs) = [10 10];

% primary input has Cload capacitance (but has no Cin)
Cload_pi = sparse(M,1);
Cload_pi(primary_inputs) = [10 10 10];

% activity frequency of gates and primary inputs
f_gates = 0.001*ones(m,1);
f_pi = sparse(M,1);
f_pi(primary_inputs) = 0.001*[10 10 10];

%********************************************************************
% derived problem data (computed from user inputs)
%********************************************************************
% fan-out cell array (compute it from the fan-in cell array)
FO = cell(M,1);
for gate = [1:m primary_outputs]
  preds = FI{gate};
  for k = 1:length(preds)
    FO{preds(k)}(end+1) = gate;
  end
end

% input and internal capacitance of gates, and driving resistance
Cin_norm  = [gates.Cin]';
Cint_norm = [gates.Cint]';
Rdrv_norm = [gates.Rdrv]';

% area specification for each gate with unit scaling
A_norm = [gates.A]';

% leakage current of gate i with unit scaling
Ileak_norm = [gates.Ileak]';

%********************************************************************
% optimization
%********************************************************************
% optimization variables
gpvar x(m)                 % scale factor
gpvar T(m)                 % arrival times

% input capacitance is an affine function of sizes
Cin  = Cin_norm.*x;
Cint = Cint_norm.*x;

% driving resistance is inversily proportional to sizes
R = Rdrv_norm./x;

% gate delay is the product of its driving resistance and load cap.
Cload = posynomial;
for gate = 1:m
  if ~ismember( FO{gate}, primary_outputs )
    Cload(gate) = sum( Cin(FO{gate}) ); 
  else
    Cload(gate) = Cin_po( FO{gate} );
  end
end
Cload = Cload';

% delay
D = 0.69*ones(m,1).*R.*( Cint + Cload );

% total area
area = A_norm'*x;

% total power calculation
Pdyn = Vdd^2*sum( f_pi(primary_inputs).*Cload_pi(primary_inputs) ) + ...
       Vdd^2*(f_gates'*(Cint + Cload));
Pstat = Vdd*Ileak_norm'*x;
power = Pdyn + Pstat;

% constraints
x_constr = ones(m,1) <= x; % all sizes greater than 1 (normalized)

% create timing constraints
timing_constr = [];
for gate = 1:m
  if ~ismember( FI{gate}, primary_inputs )
    for j = FI{gate}
      % enforce T_j + D_j <= T_i over all gates j that drive i
      timing_constr =  [timing_constr; D(gate) + T(j) <= T(gate)];
    end
  else
    % enforce D_i <= T_i for gates i connected to primary inputs
    timing_constr = [timing_constr; D(gate) <= T(gate)];
  end
end

% objective is the upper bound on the overall delay 
% and that is the max of arrival times for output gates
output_gates = [FI{primary_outputs}]; 
D = max( T(output_gates) );

% collect all the constraints
constr = [area <= Amax; power <= Pmax; timing_constr; x_constr];

% formulate the GP problem and solve it
[obj_value, solution, status] = gpsolve(D, constr);
assign(solution);

ckt_delay_plot = obj_value;
Pmax_plot = Pmax;

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
Pmax = linspace(10,20,N);
Amax = [250];
min_delay = zeros(N,1);

% set the quiet flag (no solver reporting)
global QUIET; QUIET = 1;

for n = 1:N
  % add constraints that have varying parameters
  constr(1) = area  <= Amax(k);
  constr(2) = power <= Pmax(n);

  % solve the GP problem and compute the optimal volume
  [obj_value, solution, status] = gpsolve(D, constr);
  if strcmp(status,'Infeasible')
    error('Problem is infeasible.')
  end
  min_delay(n) = obj_value;
end

% enable solver reporting again
global QUIET; QUIET = 0;

% plot the tradeoff curve

plot(Pmax,min_delay); hold on
plot(Pmax_plot, ckt_delay_plot,'o'); hold off;
xlabel('Pmax'); ylabel('Dmin');
title(['Tradeoff curve for Amax = ' num2str(Amax)])

end
% end tradeoff curve code
