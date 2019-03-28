% Runs 100 steps of fast MPC on a randomly generated system,
% using the simulation function fmpc_sim().
rand('state',0);
randn('state',0);

n = 12;
m = 3;

A = randn(n,n);
A = A./(max(abs(eig(A))));
B = rand(n,m);

% state and control limits
Xmax = 10; Umax = 0.2;
xmin = -Xmax*ones(n,1);
xmax = Xmax*ones(n,1);
umin = -Umax*ones(m,1);
umax = Umax*ones(m,1);

% objective matrices
Q = eye(n);      
R = eye(m);        

% generate process noise trajectory
nsteps = 1000;
w = 2*rand(n,nsteps)-1; w = 0.5*w;
W = (1/12)*eye(n);

% initial state
x0 = zeros(n,1);

% system description
sys.A = A ;
sys.B = B;
sys.xmax = xmax;
sys.xmin = xmin;
sys.umax = umax;
sys.umin = umin;
sys.n = n;
sys.m = m;
sys.Q = Q;
sys.R = R;

% fast MPC parameters
T = 30;                  % time horizon
params.T = T;            
params.Qf = Q;      % final state cost
params.kappa = 0.01;     % barrier parameter
params.niters = 5;       % number of newton steps
params.quiet = true;
params.nsteps = nsteps;  % number of steps to run simulation

% set up initial state and input trajectories
X = zeros(n,T);
U = zeros(m,T);
x = x0;

Qhalf = sqrtm(Q); Rhalf = sqrtm(R);
[Xhist,Uhist,telapsed] = fmpc_sim(sys,params,X,U,x,w);
Jhist = sum((Qhalf*Xhist).^2,1)+sum((Rhalf*Uhist).^2,1);
Jmean = mean(Jhist);

% report average cost and time
str = sprintf('\nThe average cost is: %f\n',Jmean);
disp(str);
str = sprintf('\nTime per step: %f\n',telapsed/nsteps);
disp(str);

% plot cost histogram
figure;
edges = linspace(0,30,50);
count = histc(Jhist,edges)';
common_axis = [0 30 0 250];

bar(edges, count, 'histc'); hold on;
hline = findobj(gca,'Type','line'); delete(hline);
hpatch = findobj(gca,'Type','patch'); 
set(hpatch,'FaceColor',[0.91,0.91,0.91]);
set(gca,'FontName','times','FontSize',16);
plot(Jmean*[1;1],[0 300],'k--');
text(Jmean+0.5,200,'Jmean');
axis([common_axis]); hold off;
title('Histogram of stage costs');


