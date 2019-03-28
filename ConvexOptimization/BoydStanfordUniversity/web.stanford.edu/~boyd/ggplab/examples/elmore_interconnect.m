% Elmore delay sizing for an interconnect network.
% (a figure is generated)
%
% This is an example taken directly from papers:
%
%   A Tutorial on Geometric Programming (see pages 31-34)
%   by Boyd, Kim, Vandenberghe, and Hassibi.
% 
%   Digital circuit optimization via geometrical programming
%   by Boyd, Kim, Patil, and Horowitz
%   Operations Research 53(6): 899-932, 2005.
%
% We consider the problem of finding optimal wire widths w_i
% of N wire segments in an interconnect network, which will
% minimize the critical Elmore delay, subject to limits on
% wire widths and the total circuit area. We use a pi-model
% for each wire segment. Problem can be formulated as GP:
%
%   minimize   D
%       s.t.   w_min <= w <= w_max
%              area  <= Amax 
%
% where variables are widths w (and arrival times T that are used
% to formulate the overall delay D expression).
%
% Important: We label root node as 1, and all the other nodes as
%            node_label_in_the_paper + 1 (due to Matlab's convention).
%            Also label nodes with increasing numbers downstream.
%
% Almir Mutapcic 02/01/2006
clear all; close all;
PLOT_TRADEOFF = 1; % to disable set PLOT_TRADEOFF = 0;

%********************************************************************
% user supplied data (problem constants and tree topology)
%********************************************************************
N = 6; % number of nodes (including the root node which is labeled as 1)

% parent node array
% specifies which node is a unique parent for node i (always have a tree)
parent(1) = 0; % root node does not have a valid parent
parent(2) = 1;
parent(3) = 2;
parent(4) = 3;
parent(5) = 2;
parent(6) = 5;

% problem constants
Rsource = 0.1;
l = 1*ones(N-1,1);
alpha = 1*ones(N-1,1);
beta  = 1*ones(N-1,1);
gamma = 1*ones(N-1,1);

% load capacitance at each node
C1 = 10; C2 = 10; C3 = 10; C4 = 10; C5 = 10;
Cload = [0 C1 C2 C3 C4 C5]; 

% minimum and maximum width and area specification 
Wmin = 1;
Wmax = 10;
Amax = 15;

%********************************************************************
% derived data (computed from user's data)
%********************************************************************
% compute children cell array (evaluate who are children for each node)
children = cell(N,1);
leafs = [];
for node = [1:N]
  children{node} = find(parent == node);
  if isempty(children{node})
    leafs(end+1) = node; % leafs have no children
  end
end

%********************************************************************
% optimization
%********************************************************************
% optimization variables
gpvar w(N-1)     % wire width
gpvar T(N)       % arrival time (Elmore delay to node i)

% wire segment resistance is inversely proportional to widths
R = alpha.*l./w;
R = [Rsource; R];

% wire segment capacitance is an affine function of widths
C_bar = beta.*l.*w + gamma.*l;
C_bar = [0; C_bar];

% compute common capacitances for each node (C_tilde in GP tutorial)
C_tilde = posynomial; % create empty posynomial
for node = [1:N]
  C_tilde(node,1) = Cload(node);
  for k = parent(node)
    if k > 0; C_tilde(node,1) = C_tilde(node,1) + C_bar(k); end;
  end
  for k = children{node}
    C_tilde(node,1) = C_tilde(node,1) + C_bar(k);
  end
end

% now compute total downstream capacitances
C_total = posynomial; % create empty posynomial
for node = N:-1:1
  C_total(node,1) = C_tilde(node);
  for k = children{node}
    C_total(node,1) = C_total(node,1) + C_total(k,1);
  end
end

% generate Elmore delay constraints
elm_delay_constr = [R(1)*C_total(1) <= T(1,1)];
for node = 2:N
  elm_delay_constr = [elm_delay_constr; ...
                      R(node)*C_total(node) + T(parent(node),1) <= T(node,1)];
end

% collect all the constraints
area = sum(w.*l);
constr(1) = area <= Amax;
constr = [constr; Wmin*ones(N-1,1) <= w; w <= Wmax*ones(N-1,1)];
constr = [constr; elm_delay_constr];

% objective is the critical Elmore delay
D = max( T(leafs) );

% solve the problem
[D_value, solution, status] = gpsolve(D, constr);
assign(solution);

% save for plotting
ckt_delay_plot = D_value;
Amax_plot = Amax;

fprintf(1,'\nOptimal Elmore delay for Amax = %d is %3.4f.\n', ...
        Amax, D_value)
disp('Optimal wire widths are: '), w

%********************************************************************
% tradeoff curve code
%********************************************************************
if( PLOT_TRADEOFF )

% set the quiet flag (no solver reporting)
global QUIET; QUIET = 1;

disp('generating the tradeoff curve')
Darray = []; Duniform = [];
for Amax = [5.05 5.25 5.5 5.75 6:25]
  % formulate the GP problem and solve it
  constr(1) = area <= Amax;
  [D_value, solution, status] = gpsolve(D, constr);
  Darray = [Darray D_value];

  % evaluate delay to leaf 4 (or 3 in the papers) with uniform sizing
  D_to_3 = C_tilde(4)*(R(1) + R(2) + R(3) + R(4)) + ...
           C_tilde(3)*(R(1) + R(2) + R(3)) + C_tilde(1)*R(1) + ...
           (C_tilde(2) + C_tilde(5) + C_tilde(6))*(R(1) + R(2));
  D_to_5 = C_tilde(6)*(R(1) + R(2) + R(5) + R(6)) + ...
           C_tilde(5)*(R(1) + R(2) + R(5)) + C_tilde(1)*R(1) + ...
           (C_tilde(2) + C_tilde(3) + C_tilde(4))*(R(1) + R(2));
  D_3_5 = max(D_to_3,D_to_5);
  Duniform = [Duniform eval(D_3_5, {'w' Amax/(N-1)*ones(N-1,1)}) ];
end

% enable solver reporting again
global QUIET; QUIET = 0;

% plot the tradeoff curve
Amax = [5.05 5.25 5.5 5.75 6:25];
plot(Darray,Amax,Duniform,Amax,'--',ckt_delay_plot,Amax_plot,'bo');
xlabel('Elmore delay D'); ylabel('Amax');
legend('optimal','uniform')

end
% end tradeoff curve code
