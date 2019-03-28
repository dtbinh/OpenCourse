% Digital circuit sizing for an inverter chain.
% (a figure is generated)
%
% This is an example taken directly from the paper:
%
%   Digital circuit optimization via geometrical programming 
%   by Boyd, Kim, Patil, and Horowitz
%   Operations Research 53(6): 899-932, 2005.
%
% We consider a chain of N inverters driving a load capacitance CL.
% The problem is to find optimal scale factors for the inverter
% that minimize the sum of them (area), while obeying constraints
% on the maximum delay through the circuit, and minimum and maximum
% limits on scale factors. There are no limits on the total power.
% (For more details about the inverter chain see sec. 2.1.11 in the paper.)
%
%   minimize   sum(x) 
%       s.t.   T_j <= Dmax          for j an output gate
%              T_j + d_i <= T_i     for j in FI(i)
%              x_min <= x <= x_max
%
% where variables are x and T.
% Here we use data structures and digital circuit models from the
% referenced paper.
%
% Almir Mutapcic 01/26/06
clear all; close all;

%********************************************************************
% problem data
%********************************************************************
N  = 8;      % number of inverters
CL = 20;     % capacitance load
Dmax = 20;   % maximum delay through the circuit
x_min = 1;   % minimum scale factor
x_max = 20;  % maximum scale factor

% primary inputs and primary outputs labels (start with N+1)
primary_inputs  = [N+1];
primary_outputs = [N+2];
M = N + length( primary_inputs ) + length( primary_outputs );

% fan-in cell array
FI{1} = [N+1];
for k = 2:N
  FI{k} = [k-1];
end
FI{N+2} = [N];

% fan-out cell array (compute it from the fan-in cell array)
FO = cell(M,1);
for gate = [1:N primary_outputs]
  preds = FI{gate};
  for k = 1:length(preds)
    FO{preds(k)}(end+1) = gate;
  end
end

% input and internal capacitance of gates and the driving resistance
Cin_norm  = ones(N,1);
Cint_norm = ones(N,1);
Rdrv_norm = ones(N,1);

% place extra capacitance before the input of the 5th inverter
Cin_norm(5) = 80;

% primary output has Cin capacitance (but has no Cload)
Cin_po = sparse(M,1);
Cin_po(primary_outputs) = CL;

% primary input has Cload capacitance (but has no Cin)
Cload_pi = sparse(M,1);
Cload_pi(primary_inputs) = 1;

%********************************************************************
% optimization
%********************************************************************
% optimization variables
gpvar x(N)                 % sizes
gpvar T(N)                 % arrival times

% input capacitance is an affine function of sizes
Cin  = Cin_norm.*x;
Cint = Cint_norm.*x;

% driving resistance is inversily proportional to sizes
R = Rdrv_norm./x;

% gate delay is the product of its driving resistance and load cap.
Cload = posynomial;
for gate = 1:N
  if ~ismember( FO{gate}, primary_outputs )
    Cload(gate) = sum( Cin(FO{gate}) ); 
  else
    Cload(gate) = Cin_po( FO{gate} );
  end
end
Cload = Cload';

% delay
D = 0.69*ones(N,1).*R.*( Cint + Cload );

% create timing constraints
timing_constr = [];
for gate = 1:N
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

% circuit delay is the max of arrival times for output gates
output_gates = [FI{primary_outputs}]; 
circuit_delay = max( T(output_gates) );

% collect all the constraints
constr = [timing_constr; circuit_delay <= Dmax; 
          x_min*ones(N,1) <= x; x <= x_max*ones(N,1)];

% minimize the sum of scale factors subject to above constraints
[obj_value, solution, status] = gpsolve(sum(x), constr);
assign(solution);

% message about extra capacitance and result display
disp(' ')
disp(['Note: there is an extra capacitance between the 4th and 5th inverter'...
     ' in the chain.'])
fprintf(1,'\nOptimal scale factors are: \n'), x

% plot scale factors and maximum delay for inverter i
close all;
subplot(2,1,1); plot([1:N],T,'g--',[1:N],T,'bo');
ylabel('maximum delay T')
subplot(2,1,2); stem([1:N],x);
ylabel('scale factor x')
xlabel('inverter stage')
