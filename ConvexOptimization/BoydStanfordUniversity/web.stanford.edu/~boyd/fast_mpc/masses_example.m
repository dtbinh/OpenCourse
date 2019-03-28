% Runs 100 steps of fast MPC on an example with 6 masses
% connected by springs, using fmpc_step().
rand('state',0); % make
k = 1;          % spring constant 
lam = 0;        % damping constant 
a = -2*k;
b = -2*lam;
c = k;
d = lam;

n = 12; % state dimension
m = 3; % input dimension

Acts = [zeros(n/2) eye(n/2);
    [a,c,0,0,0,0,b,d,0,0,0,0;
     c,a,c,0,0,0,d,b,d,0,0,0;
     0,c,a,c,0,0,0,d,b,d,0,0;
     0,0,c,a,c,0,0,0,d,b,d,0;
     0,0,0,c,a,c,0,0,0,d,b,d;
     0,0,0,0,c,a,0,0,0,0,d,b]];

Bcts = [zeros(n/2,m);
    [1, 0, 0;
    -1, 0, 0;
     0, 1, 0;
     0, 0, 1;
     0,-1, 0;
     0, 0,-1]];

% convert to discrete-time system
ts = 0.5;       % sampling time
A = expm(ts*Acts);
B = (Acts\(expm(ts*Acts)-eye(n)))*Bcts;

% objective matrices
Q = eye(n);      
R = eye(m);        

% state and control limits
Xmax = 4; Umax = 0.5;
xmin = -Xmax*ones(n,1);
xmax = Xmax*ones(n,1);
umin = -Umax*ones(m,1);
umax = Umax*ones(m,1);

% process noise trajectory
nsteps = 100;
w = 2*rand(n,nsteps)-1; w(1:n/2,:) = 0; w = 0.5*w;
W = (1/12)*diag([0;0;0;0;0;0;1;1;1;1;1;1]);

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
T = 10;                % time horizon
params.T = T; 
params.Qf = Q;         % final state cost
params.kappa = 0.01;   % barrier parameter
params.niters = 5;     % number of newton steps
params.quiet = false;

% allocate history matrices
Xhist = zeros(n,nsteps);  % state
Uhist = zeros(m,nsteps);  % input
Jhist = zeros(1,nsteps);  % stage cost
thist = zeros(1,nsteps);  % fmpc run time

% set up initial state and input trajectories
X = zeros(n,T);
U = zeros(m,T);
x = x0;

tic;
for i = 1:nsteps
    [X,U,telapsed] = fmpc_step(sys,params,X,U,x);
    u = U(:,1);
    % record state, input, stage cost, and fmpc run time
    Xhist(:,i) = x; Uhist(:,i) = u;
    Jhist(i) = x'*Q*x+u'*R*u;
    thist(i) = telapsed;
    % state update
    x = A*x+B*u+w(:,i);
    % shift previous state and input trajectories 
    % for warm start in next step
    X = [X(:,2:T),zeros(n,1)];
    U = [U(:,2:T),zeros(m,1)];
end
t = toc;

% print average cost 
str = sprintf('\nThe average cost is: %f\n',mean(Jhist));
disp(str);

% plot trajectories for x_1 and u_1
tvec = 0:0.5:49.5; figure;
subplot(2,1,1); stairs(tvec,Xhist(1,1:100));
axis([0,50,-2,2]); ylabel('x_1');
subplot(2,1,2); stairs(tvec,Uhist(1,1:100));
axis([0,50,-0.6,0.6]);
xlabel('t'); ylabel('u_1');
