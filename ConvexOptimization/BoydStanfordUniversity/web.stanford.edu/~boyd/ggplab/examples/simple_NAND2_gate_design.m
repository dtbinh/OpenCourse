% Two-input NAND gate sizing (a simple, hard-coded example).
% (a figure is generated if the tradeoff flag is turned on)
%
% This is an example taken directly from the paper:
%
%   Digital circuit optimization via geometrical programming
%   by Boyd, Kim, Patil, and Horowitz
%   Operations Research 53(6): 899-932, 2005.
%
% Solves the problem of choosing device widths w_i for the given
% NAND2 gate in order to achive minimum Elmore delay for different
% gate transitions, subject to limits on the device widths, 
% gate area, power, and so on. The problem is a GP:
%
%   minimize   D = max( D_1, ..., D_k )  for k transitions
%       s.t.   w_min <= w <= w_max
%              A <= Amax, etc.
%
% where variables are widths w.
%
% This code is specific to the NAND2 gate shown in figure 19 
% (page 926) of the paper. All the constraints and the objective
% are hard-coded for this particular circuit.
% 
% Almir Mutapcic 02/01/2006
clear all; close all;
PLOT_TRADEOFF = 1; % to disable set PLOT_TRADEOFF = 0;

%********************************************************************
% problem data and hard-coded GP specification (evaluate all transitions)
%********************************************************************
N = 4;       % number of devices
Cload = 12;  % load capacitance
Vdd = 1.5;   % voltage

% device specs
NMOS = struct('R',0.4831, 'Cdb',0.6, 'Csb',0.6, 'Cgb',1, 'Cgs',1);
PMOS = struct('R',2*0.4831, 'Cdb',0.6, 'Csb',0.6, 'Cgb',1, 'Cgs',1);

% device width variables
gpvar w(N)

% gate specs
gates(1:2) = PMOS; gates(3:4) = NMOS;

for num = 1:N
  gates(num).R   = gates(num).R/w(num);
  gates(num).Cdb = gates(num).Cdb*w(num);
  gates(num).Csb = gates(num).Csb*w(num);
  gates(num).Cgb = gates(num).Cgb*w(num);
  gates(num).Cgs = gates(num).Cgs*w(num);
end

% capacitances
C1 = sum([gates(1:3).Cdb]) + Cload;
C2 = gates(3).Csb + gates(4).Cdb;

% input capacitances
Cin_A = sum([ gates([2 3]).Cgb ]) + sum([ gates([2 3]).Cgs ]);
Cin_B = sum([ gates([1 4]).Cgb ]) + sum([ gates([1 4]).Cgs ]);

% resistances
R = [gates.R]';

% delays and dissipated energies for all six possible transitions 
% transition 1 is A: 1->1, B: 1->0, Z: 0->1
D1 = R(1)*(C1 + C2);
E1 = (C1 + C2)*Vdd^2/2;
% transition 2 is A: 1->0, B: 1->1, Z: 0->1
D2 = R(2)*C1;
E2 = C1*Vdd^2/2;
% transition 3 is A: 1->0, B: 1->0, Z: 0->1
% D3 = C1*R(1)*R(2)/(R(1) + R(2)); % not a posynomial
E3 = C1*Vdd^2/2;
% transition 4 is A: 1->1, B: 0->1, Z: 1->0
D4 = C1*R(3) + R(4)*(C1 + C2);
E4 = (C1 + C2)*Vdd^2/2;
% transition 5 is A: 0->1, B: 1->1, Z: 1->0
D5 = C1*(R(3) + R(4));
E5 = (C1 + C2)*Vdd^2/2;
% transition 6 is A: 0->1, B: 0->1, Z: 1->0
D6 = C1*R(3) + R(4)*(C1 + C2);
E6 = (C1 + C2)*Vdd^2/2;

% maximum area and power specification 
Amax = [5:45];
Amax = 24;
wmin = 1;
area = sum(w);

% objective is the worst-case delay
gate_delay = max(D1,D2,D4); 

% collect all the constraints
constr = [area <= Amax; wmin*ones(N,1) <= w];

% formulate the GP problem and solve it
[obj_value, solution, status] = gpsolve(gate_delay, constr);
assign(solution);

fprintf(1,'\nOptimal gate delay for Amax = %d is %3.4f.\n', ...
        Amax, obj_value)
disp('Optimal device widths are: '), w

%********************************************************************
% tradeoff curve code
%********************************************************************
if( PLOT_TRADEOFF )

% varying parameters for an optimal trade-off curve
Npoints = 25;
Amax = linspace(5,45,Npoints);
Dopt = []; Duniform = [];

% set the quiet flag (no solver reporting)
global QUIET; QUIET = 1;

for k = 1:Npoints
  % add constraints that have varying parameters
  constr(1) = area <= Amax(k);

  % solve the GP problem and compute the optimal volume
  [obj_value, solution, status] = gpsolve(gate_delay, constr);
  Dopt = [Dopt obj_value];
  Duniform = [Duniform eval(gate_delay, {'w' Amax(k)/N*ones(N,1)}) ];
end

% enable solver reporting again
global QUIET; QUIET = 0;

% plot the tradeoff curve
plot(Dopt,Amax,Duniform,Amax,'--');
xlabel('Dmin'); ylabel('Amax');
legend('optimal','uniform')

end
% end tradeoff curve code
