% Matlab script for the paper
% "Optimized Slowdown in Real-Time Task Systems via Geometric Programming",
% by A. Mutapcic, S. Murali, S. Boyd, R. Gupta, D. Atienza, G. De Micheli.
%
% This script requires gpposy, a Matlab interior-point GP solver,
% which is available from the same website as this script, or
% as part of the ggplab toolbox: http://www.stanford.edu/~boyd/ggplab/

%********************************************************************
% create problem data (or load your own T, D, C task timing info)
%********************************************************************
fprintf(1,'Loading problem constants and task time specifications.\n');
fprintf(1,'We have:\n');

% number of tasks
n = 20;

% randomly generate task timing info (period, deadline, WCET)
rand('state',1);
T = ceil(random('unif',20,50,n,1))*1000;   % period
D = ceil(T/4);                             % deadline
C = ceil(random('unif',100,500,n,1)*10/n); % WCET (worst-case exec time)

% problem constants
Vth = 0.36;      % treshold voltage
Vddmax = 1.8;    % maximum operating voltage is 1.8 V at 1 GHz
alpha = 1.5;
k2 = (Vddmax-Vth)^alpha/(Vddmax);
N = ones(n,1);

% (quasi) hyper-period approximation
epsilon = 0.01;
Tmax = (1/epsilon)*(max(T-D));

fprintf(1,'   n    is %g\n', n);
fprintf(1,'   Tmax is %g\n', Tmax);

fprintf(1,'Processing timing constraints: ');
S = [];
for i=1:n
  S = union(S,[D(i):T(i):Tmax]);
end
L = length(S);

% compute param
param = zeros(n,L);
for i=1:L
  param(:,i) = C/S(i).*(floor((S(i) - D)./T) + 1);
end
fprintf(1,'we have %d constraints.\n',L); 

%********************************************************************
% custom parsing for gpposy
%********************************************************************
fprintf(1,'Parsing gpposy data, constraint number: \n  '); 
I = speye(n); Z = sparse(n,n);

% objective function
A0 = [2*I, Z]; b0 = N; szs = [n];

% timing constraints
A1 = []; b1 = [];
for t = 1:size(param,2)
  % equation: param(:,t)'*(eta.^-1) <= Tnew(t);
  A1 = [A1; Z, -I]; b1 = [b1; param(:,t)]; szs = [szs; n];
  if rem(t,100) == 0, fprintf(1,' %d',t), end
end
fprintf(1,'\n');

% voltage-slowdown relation
A2 = []; b2 = [];
for k = 1:n
  % equation: k2*(Vth+dV).*dV.^(-alpha) <= eta.^-1;
  mtx = sparse(2,2*n);
  mtx(1,k) = 1/alpha - 1; mtx(1,k+n) = 1/alpha;
  mtx(2,k) = -1;
  A2 = [A2; mtx]; b2 = [b2; k2^(1/alpha); Vth]; szs = [szs; 2];
end

% lower and upper bounds
l = [Vth*ones(n,1); zeros(n,1)];
u = [Vddmax*ones(n,1); ones(n,1)];

fprintf(1,'Processed GP problem: total matrix row length is %d.\n',sum(szs));

%********************************************************************
% call gpposy (GP solver)
%********************************************************************
EPS = 10^-8;
quiet = 0; % quiet = 1;
A = [A0; A1; A2]; b = [b0; b1+EPS; b2];
fprintf(1,'Calling gpposy to solve the problem.\n');
tic
[x,status,lambda,nu] = gpposy(A,b,szs,[],[],l,u,quiet);
toc

if strcmp( status, 'Solved') 
  V_opt = x(1:n);
  eta_opt = x(n+1:2*n);
  fprintf(1,'\nOptimal objective value (min energy) is %3.4f.\n',N'*V_opt.^2);
  fprintf(1,'Minimum slowdown factor is %3.4f.\n',min(eta_opt));
  fprintf(1,'To display optimal variables, type >> V_opt, eta_opt\n\n');
else
  fprintf(1,'**Problem is %s**.\n',status);
end

%
% the following code computes the optimized slowdown using CVX toolbox
% (for comparison and verification purposes)
%
% Warning: only enable this code for problems with very small n (2-5),
% and small epsilon (0.1), otherwise it takes a very long time to parse
% the problem and start the computation
%
if( 0 ) 
fprintf(1,'\n\nUsing CVX to verify results.\n');
%********************************************************************
% CVX code
%********************************************************************
cvx_begin gp
variables dV(n) eta(n)

minimize( N'*((Vth+dV).^2) )
subject to
  
  for t = 1:L
    t
    param(:,t)'*(eta.^-1) <= 1;
  end
    
  k2*(Vth+dV).*dV.^(-alpha) <= eta.^-1;
  eta <= 1;
cvx_end

% result comparisons
fprintf(1,'Result comparisons:\n')
V_opt, V_cvx = Vth + dV
eta_opt, eta_cvx = eta

end % CVX
