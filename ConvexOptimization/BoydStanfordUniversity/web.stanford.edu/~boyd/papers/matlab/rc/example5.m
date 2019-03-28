%
% EXAMPLE 5:  Combined wires sizing and spacing
%

% circuit parameters
n = 6;         % number of nodes per wire
N = 3*n;       % total number of nodes
m = n-1;       % number of segments per wire 
alpha = 1;     % conductance per segment is is alpha*size 
beta = 0.5;    % capacitance per segment is twice beta*size 
gamma = 2;     % coupling capacitance is twice gamma*distance
G0 = 100;      % source conductance
C01 = 10;      % load of first wire
C02 = 20;      % load of second wire
C03 = 30;      % load of third wire    
wmin = 0.1;    % minimum width
wmax = 2.0;    % maximum width 
smin = 1.0;    % minimum distance between wires
smax = 50;     % upper bound on s1 and s2  (meant to be inactive)


% variables:
% x(1:m):      sizes of first wire  (w1 in paper)
% x(m+1:2m):   sizes of second wire (w2 in paper)
% x(2m+1:3m):  sizes of third wire (w3 in paper)
% x(3m+1:4m):  inverse of distances between wire 1 and 2 (t1 in paper)
% x(4m+1:5m):  inverse of distances between wire 2 and 3 (t2 in paper)
% x(5m+1):     s1
% x(5m+2):     s2


% capacitance matrix

CC = zeros(N*N, 5*m+3);

% constant terms
CC(N*(n-1)+n,1) = C01;
CC(N*(2*n-1)+2*n,1) = C02;
CC(N*(3*n-1)+3*n,1) = C03;

% capacitances to ground
for i=1:m 

   % first wire
   CC([(i-1)*N+i, i*N+i+1], i+1) = beta*[1;1];

   % second wire
   CC([(n+i-1)*N+n+i, (n+i)*N+n+i+1], m+i+1) = beta*[1;1];

   % third wire
   CC([(2*n+i-1)*N+2*n+i, (2*n+i)*N+2*n+i+1], 2*m+i+1) = beta*[1;1];

end;

% coupling capacitors

for i=1:m

   % between first and second wire
   CC([(i-1)*N+i, (i-1)*N+n+i, ...
       (n+i-1)*N+i, (n+i-1)*N+n+i], 3*m+1+i) ...
   = CC([(i-1)*N+i, (i-1)*N+n+i, ...
         (n+i-1)*N+i, (n+i-1)*N+n+i], 3*m+1+i) + gamma*[1;-1;-1;1];

   CC([i*N+i+1, i*N+n+i+1, (n+i)*N+i+1, (n+i)*N+n+i+1], 3*m+1+i) ...
   = CC([i*N+i+1, i*N+n+i+1, (n+i)*N+i+1, (n+i)*N+n+i+1], 3*m+1+i) ...
      + gamma*[1;-1;-1;1];

   % between second and third wire
   CC([(n+i-1)*N+n+i, (n+i-1)*N+2*n+i, (2*n+i-1)*N+n+i, ...
       (2*n+i-1)*N+2*n+i], 4*m+1+i) ...
   = CC([(n+i-1)*N+n+i, (n+i-1)*N+2*n+i, (2*n+i-1)*N+n+i, ...
         (2*n+i-1)*N+2*n+i], 4*m+1+i) + gamma*[1;-1;-1;1];

   CC([(n+i)*N+n+i+1, (n+i)*N+2*n+i+1, ...
       (2*n+i)*N+n+i+1, (2*n+i)*N+2*n+i+1], 4*m+1+i) ...
   = CC([(n+i)*N+n+i+1, (n+i)*N+2*n+i+1, ...
         (2*n+i)*N+n+i+1, (2*n+i)*N+2*n+i+1], 4*m+1+i) ...
     + gamma*[1;-1;-1;1];

end;


% conductance matrix

GG = zeros(N*N, 5*m+3);

% constant terms
GG(1,1) = G0; 
GG(n*N+n+1,1) = G0;
GG(2*n*N+2*n+1,1) = G0;

% segment conductances
for i=1:m

   % first wire
   GG([(i-1)*N+i, (i-1)*N+i+1, i*N+i, i*N+i+1], i+1) ...
   = GG([(i-1)*N+i, (i-1)*N+i+1, i*N+i, i*N+i+1], i+1) ...
     + alpha*[1;-1;-1;1];

   % second wire
   GG([(n+i-1)*N+n+i, (n+i-1)*N+n+i+1, ...
       (n+i)*N+n+i, (n+i)*N+n+i+1], m+i+1) ...
   = GG([(n+i-1)*N+n+i, (n+i-1)*N+n+i+1, ...
         (n+i)*N+n+i, (n+i)*N+n+i+1], m+i+1) + alpha*[1;-1;-1;1];

   % third wire
   GG([(2*n+i-1)*N+2*n+i, (2*n+i-1)*N+2*n+i+1, ...
       (2*n+i)*N+2*n+i, (2*n+i)*N+2*n+i+1], 2*m+i+1) ...
   = GG([(2*n+i-1)*N+2*n+i, (2*n+i-1)*N+2*n+i+1, ...
         (2*n+i)*N+2*n+i, (2*n+i)*N+2*n+i+1], 2*m+i+1) ...
     + alpha*[1;-1;-1;1];

end;


% LMI constraints
% 
%  G(x) - (1/delay)*C(x) (+t*I) >= 0   (LMI of size NxN)
%  
%  [(t1i+s1-w1i-0.5*w2i) 0                    2                    ]
%  [0                    (t1i+s1-w1i-0.5*w2i) (t1i-s1+w1i+0.5*w2i) ]
%  [2                    (t1i-s1+w1i+0.5*w2i) (t1i+s1-w1i-0.5*w2i) ]
%     (+t*I) >= 0 for i=1,...,m        (m LMIs of size 3x3) 
%
%  [(t2i+s2-w3i-0.5*w2i) 0                    2                    ]
%  [0                    (t2i+s2-w3i-0.5*w2i) (t2i-s2+w3i+0.5*w2i) ]
%  [2                    (t2i-s2+w3i+0.5*w2i) (t2i+s2-w3i-0.5*w2i) ]
%     (+t*I) >= 0 for i=1,...,m        (2m LMIs of size 3x3) 
%
%  s1 - w1 - 0.5*w2 (+t) .> smin; 
%  s2 - w3 - 0.5*w2 (+t) .> smin;      (2m scalars) 
%
%  w1,w2,w3 >= wmin
%  t1,t2  >= 0
%  s1,s2  >= 0                         (5m+2 scalars)
%
%  w1,w2,w3 <= wmax
%  t1,t2  <= 1/smin 
%  s1,s2  <= smax                      (5m+2 scalars)
%
 
% 5m+2 variables:
%
% x(1:m):      sizes of first wire  (w1 in paper)
% x(m+1:2m):   sizes of second wire (w2 in paper)
% x(2m+1:3m):  sizes of third wire (w3 in paper)
% x(3m+1:4m):  inverse of distances between wire 1 and 2 (t1 in paper)
% x(4m+1:5m):  inverse of distances between wire 2 and 3 (t2 in paper)
% x(5m+1):     s1
% x(5m+2):     s2


FF = zeros(N*N + 2*m*9 + 2*m + 2*(5*m+2), 5*m+4);
blck_szs = [N; 3*ones(2*m,1); ones(12*m+4,1)]; 

% identity matrix in last column
FF(1:N*N,5*m+4) = reshape(eye(N),N*N,1);
for i=1:m
   FF(N*N+(i-1)*9+[1:9],5*m+4) = reshape(eye(3),9,1);
   FF(N*N+m*9+(i-1)*9+[1:9],5*m+4) = reshape(eye(3),9,1);
end
FF(N*N+m*18+[1:2*m],5*m+4) = ones(2*m,1);

% hyperbolic constraints
%  [(t1i+s1-w1i-0.5*w2i) 0                    2                    ]
%  [0                    (t1i+s1-w1i-0.5*w2i) (t1i-s1+w1i+0.5*w2i) ]
%  [2                    (t1i-s1+w1i+0.5*w2i) (t1i+s1-w1i-0.5*w2i) ]
%     (+t*I) >= 0 for i=1,...,m        (m LMIs of size 3x3) 
%
%  [(t2i+s2-w3i-0.5*w2i) 0                    2                    ]
%  [0                    (t2i+s2-w3i-0.5*w2i) (t2i-s2+w3i+0.5*w2i) ]
%  [2                    (t2i-s2+w3i+0.5*w2i) (t2i+s2-w3i-0.5*w2i) ]
%     (+t*I) >= 0 for i=1,...,m        (2m LMIs of size 3x3) 

for i=1:m

   % constant term in position 3 and 7
   FF(N*N+(i-1)*9+[3,7],1) = [2;2];
   FF(N*N+m*9+(i-1)*9+[3,7],1) = [2;2];

   % s1, s2 term
   FF(N*N+(i-1)*9+[1,5,6,8,9],5*m+2) = [1;1;-1;-1;1];
   FF(N*N+m*9+(i-1)*9+[1,5,6,8,9],5*m+3) = [1;1;-1;-1;1];

   % w1 terms
   FF(N*N+(i-1)*9+[1,5,6,8,9],1+i) = [-1;-1;1;1;-1];

   % w2 terms
   FF(N*N+(i-1)*9+[1,5,6,8,9],1+m+i) = 0.5*[-1;-1;1;1;-1];
   FF(N*N+m*9+(i-1)*9+[1,5,6,8,9],1+m+i) = 0.5*[-1;-1;1;1;-1];

   % w3 terms
   FF(N*N+m*9+(i-1)*9+[1,5,6,8,9],1+2*m+i) = [-1;-1;1;1;-1];

   % t1 terms
   FF(N*N+(i-1)*9+[1,5,6,8,9],1+3*m+i) = [1;1;1;1;1];

   % t2 terms
   FF(N*N+m*9+(i-1)*9+[1,5,6,8,9],1+4*m+i) = [1;1;1;1;1];

end;

%  s1 - w1 - 0.5*w2 .> smin; 
%  s2 - w3 - 0.5*w2 .> smin;    (2m scalars) 

FF(N*N+m*18+[1:m],1+[1:m]) = -eye(m);          %w1
FF(N*N+m*18+[1:m],m+1+[1:m]) = -0.5*eye(m);    %w2
FF(N*N+m*18+[1:m],5*m+2) = ones(m,1);          %s1

FF(N*N+m*18+m+[1:m],1+2*m+[1:m]) = -eye(m);    %w3
FF(N*N+m*18+m+[1:m],m+1+[1:m]) = -0.5*eye(m);  %w2
FF(N*N+m*18+m+[1:m],5*m+3) = ones(m,1);        %s1

% lower bounds 
FF(N*N+m*18+2*m+[1:5*m+2],1+[1:5*m+2]) = eye(5*m+2);
FF(N*N+m*18+2*m+5*m+2+[1:5*m+2],1+[1:5*m+2]) = -eye(5*m+2);

% cost function
c = [zeros(5*m,1); 1; 1];



%%%%%%%%%%  compute tradeoff curve  %%%%%%%%%%

% tolerances for SP
abstol = 1e-4;
reltol = 1e-4;

% points on tradeoff curve
nopts = 10;
delays = linspace(85,200,nopts);
areas = zeros(1,nopts);

% compute tradeoff curve
for j=1:nopts
   
   delay = delays(j);
   disp([' ']);
   disp(['Point ', int2str(j), ...
        ' on the tradeoff curve (delay = ', num2str(delay), ').']);

   % first block of F(x) = G(x) - (1/delay)*C(x)
   FF(1:N*N,1:5*m+3) = GG - (1/delay)*CC;

   % try very thick wire
   x0 = 0.99*[wmax*ones(3*m,1); (1/smin)*ones(2*m,1); smax; smax];  
   t = min(eig(reshape(FF(1:N*N,1:5*m+3)*[1;x0],N,N)));
   for i=1:2*n
      t = min(t, min(eig(reshape(FF(N*N+(i-1)*9+[1:9],1:5*m+3) ...
         * [1;x0],3,3))));
   end;

   feasible = 1;
   if (t <= 0)  % x0 is infeasible, run phase 1

      disp([' ']);  disp(['Phase 1.']);

      % use slightly tighter bounds on x to ensure strict feasibility
      FF(N*N+18*m+2*m+[1:3*m],1) = -(1-1e-5)*wmin*ones(3*m,1);
      FF(N*N+18*m+2*m+3*m+[1:2*m+2],1) = -1e-5*ones(2*m+2,1);

      FF(N*N+18*m+2*m+5*m+2+[1:3*m],1) = (1-1e-5)*wmax*ones(3*m,1);
      FF(N*N+18*m+2*m+8*m+2+[1:2*m],1) = (1-1e-5)*(1/smin)*ones(2*m,1);
      FF(N*N+18*m+2*m+10*m+2+[1:2],1) = (1-1e-5)*smax*ones(2,1);
     
      % dual feasible solution
      Z0 = zeros(sum(blck_szs.^2),1);
      Z0(1:N*N) = reshape(eye(N)/(N+6*m+2*m),N*N,1);
      for i=1:2*m
         Z0(N*N+9*(i-1)+[1:9]) = reshape(eye(3)/(N+6*m+2*m),9,1);
      end;
      Z0(N*N+18*m+[1:2*m]) = ones(2*m,1)/(N+6*m+2*m);
      z = FF(1:(N*N+18*m+2*m),2:5*m+3)' * Z0(1:N*N+18*m+2*m);
      Z0(N*N+18*m+2*m+[1:10*m+4]) = [-min(z,0)+1e-3; max(z,0)+1e-3 ];

      % call SP
      [x,Z,ul,info,time] = sp(FF, blck_szs, [zeros(5*m+2,1);1], ...
          [x0; -1.1*t], Z0, 10.0, abstol, -1.0, 0.0, 100);

      if (ul(1) >= 0),
         areas(i) = Inf;
         if (ul(2) >= 0),
            disp(['Infeasible.']);
         else
            disp(['Feasibility could not be determined.']);
         end;
         feasible = 0;
      else x0 = x(1:5*m+2);    % initial point for phase 2.
      end;

   end;

   if (feasible) % use x0 as starting point, run phase 2

      disp([' ']); disp(['Phase 2.']);

      % use original bounds on x again
      FF(N*N+18*m+2*m+[1:3*m],1) = -wmin*ones(3*m,1);
      FF(N*N+18*m+2*m+3*m+[1:2*m+2],1) = zeros(2*m+2,1);

      FF(N*N+18*m+2*m+5*m+2+[1:3*m],1) = wmax*ones(3*m,1);
      FF(N*N+18*m+2*m+8*m+2+[1:2*m],1) = (1/smin)*ones(2*m,1);
      FF(N*N+18*m+2*m+10*m+2+[1:2],1) = smax*ones(2,1);

      % dual feasible point
      Z0 = zeros(sum(blck_szs.^2),1);
      Z0(1:N*N) = reshape(eye(N)/(N+6*m+2*m),N*N,1);
      for i=1:2*m
         Z0(N*N+9*(i-1)+[1:9]) = reshape(eye(3)/(N+6*m+2*m),9,1);
      end;
      Z0(N*N+18*m+[1:2*m]) = ones(2*m,1)/(N+6*m+2*m);
      z = FF(1:(N*N+18*m+2*m),2:5*m+3)' * Z0(1:N*N+18*m+2*m) - c;
      Z0(N*N+18*m+2*m+[1:10*m+4]) = [-min(z,0)+1e-3; max(z,0)+1e-3 ];

      % call SP
      [x,Z,ul,info,time] = sp(FF(:,1:5*m+3), blck_szs, c, ...
          x0, Z0, 10.0, abstol, reltol, 0.0, 100);

      % optimal value
      areas(j) = c'*x;

   end;

end;


figure(1);
ind = finite(areas);
plot(areas(ind), delays(ind));
xlabel('area');
ylabel('Tdom');



%%%%%%%%%%  compute two solutions  %%%%%%%%%%

% tolerances for SP
abstol = 1e-6;
reltol = 1e-6;

% two points on the tradeoff curve
delays = [127 87];
sizes = zeros(5*m+2,2);

for j=[1,2]
   
   delay = delays(j);
   disp([' ']);
   disp(['Compute solution ', int2str(j), ...
        ' (delay = ', num2str(delay), ').']);

   % first block of F(x) = G(x) - (1/delay)*C(x)
   FF(1:N*N,1:5*m+3) = GG - (1/delay)*CC;

   % try very thick wire
   x0 = 0.99*[wmax*ones(3*m,1); (1/smin)*ones(2*m,1); smax; smax];  
   t = min(eig(reshape(FF(1:N*N,1:5*m+3)*[1;x0],N,N)));
   for i=1:2*n
      t = min(t, min(eig(reshape(FF(N*N+(i-1)*9+[1:9],1:5*m+3) ...
         * [1;x0],3,3))));
   end;

   feasible = 1;
   if (t <= 0)  % x0 is infeasible, run phase 1

      disp([' ']);  disp(['Phase 1.']);

      % use slightly tighter bounds on x to ensure strict feasibility
      FF(N*N+18*m+2*m+[1:3*m],1) = -(1-1e-5)*wmin*ones(3*m,1);
      FF(N*N+18*m+2*m+3*m+[1:2*m+2],1) = -1e-5*ones(2*m+2,1);

      FF(N*N+18*m+2*m+5*m+2+[1:3*m],1) = (1-1e-5)*wmax*ones(3*m,1);
      FF(N*N+18*m+2*m+8*m+2+[1:2*m],1) = (1-1e-5)*(1/smin)*ones(2*m,1);
      FF(N*N+18*m+2*m+10*m+2+[1:2],1) = (1-1e-5)*smax*ones(2,1);
     
      % dual feasible solution
      Z0 = zeros(sum(blck_szs.^2),1);
      Z0(1:N*N) = reshape(eye(N)/(N+6*m+2*m),N*N,1);
      for i=1:2*m
         Z0(N*N+9*(i-1)+[1:9]) = reshape(eye(3)/(N+6*m+2*m),9,1);
      end;
      Z0(N*N+18*m+[1:2*m]) = ones(2*m,1)/(N+6*m+2*m);
      z = FF(1:(N*N+18*m+2*m),2:5*m+3)' * Z0(1:N*N+18*m+2*m);
      Z0(N*N+18*m+2*m+[1:10*m+4]) = [-min(z,0)+1e-3; max(z,0)+1e-3 ];

      % call SP
      [x,Z,ul,info,time] = sp(FF, blck_szs, [zeros(5*m+2,1);1], ...
          [x0; -1.1*t], Z0, 10.0, abstol, -1.0, 0.0, 100);

      if (ul(1) >= 0),
         if (ul(2) >= 0),
            disp(['Infeasible.']);
         else
            disp(['Feasibility could not be determined.']);
         end;
         feasible = 0;
      else x0 = x(1:5*m+2);    % initial point for phase 2.
      end;

   end;

   if (feasible) % use x0 as starting point, run phase 2

      disp([' ']); disp(['Phase 2.']);

      % use original bounds on x again
      FF(N*N+18*m+2*m+[1:3*m],1) = -wmin*ones(3*m,1);
      FF(N*N+18*m+2*m+3*m+[1:2*m+2],1) = zeros(2*m+2,1);

      FF(N*N+18*m+2*m+5*m+2+[1:3*m],1) = wmax*ones(3*m,1);
      FF(N*N+18*m+2*m+8*m+2+[1:2*m],1) = (1/smin)*ones(2*m,1);
      FF(N*N+18*m+2*m+10*m+2+[1:2],1) = smax*ones(2,1);

      % dual feasible point
      Z0 = zeros(sum(blck_szs.^2),1);
      Z0(1:N*N) = reshape(eye(N)/(N+6*m+2*m),N*N,1);
      for i=1:2*m
         Z0(N*N+9*(i-1)+[1:9]) = reshape(eye(3)/(N+6*m+2*m),9,1);
      end;
      Z0(N*N+18*m+[1:2*m]) = ones(2*m,1)/(N+6*m+2*m);
      z = FF(1:(N*N+18*m+2*m),2:5*m+3)' * Z0(1:N*N+18*m+2*m) - c;
      Z0(N*N+18*m+2*m+[1:10*m+4]) = [-min(z,0)+1e-3; max(z,0)+1e-3 ];

      % call SP
      [sizes(:,j),Z,ul,info,time] = sp(FF(:,1:5*m+3), blck_szs, c, ...
          x0, Z0, 10.0, abstol, reltol, 0.0, 100);

      % optimal value

   end;

end;


%%%%%%%%%%  show first solution  %%%%%%%%%%

figure(2); 
hold off;
colormap(gray);

w1 = sizes(1:m,1);
w2 = sizes(m+[1:m],1);
w3 = sizes(2*m+[1:m],1);
s1 = sizes(5*m+1,1);
s2 = sizes(5*m+2,1);

x = zeros(10,1);
x([1:2:2*m-1],:) = [0:m-1]';
x([2:2:2*m],:) = [1:m]';

width1 = zeros(2*m,1);
width1([1:2:2*m-1]) = s1+s2-w1;
width1([2:2:2*m]) = s1+s2-w1;

width2a = zeros(2*m,1);
width2a([1:2:2*m-1]) = s2+w2/2;
width2a([2:2:2*m]) = s2+w2/2;
width2b = zeros(2*m,1);
width2b([1:2:2*m-1]) = s2-w2/2;
width2b([2:2:2*m]) = s2-w2/2;

width3 = zeros(2*m,1);
width3([1:2:9]) = w3;
width3([2:2:10]) = w3;

plot([x;flipud(x);0], [width1;(s1+s2)*ones(size(x));width1(1)],'-', ...
     [x;flipud(x);0], [width2a;flipud(width2b);width2a(1)], '-', ...
     [x;flipud(x);0], [width3;zeros(size(x));width3(1)], '-');
axis([-0.1, m+0.1,-0.1, s1+s2+0.1]);
hold on
fill([x;m;0]',[width1;s1+s2;s1+s2]', 0.9*ones(size([x;m;0]')));
fill([x;flipud(x)]',[width2a;flipud(width2b)]', ...
     0.9*ones(size([x;x]')));
fill([x;m;0]',[width3;0;0]', 0.9*ones(size([x;0;0]')));
caxis([-1,1])
plot([x;flipud(x);0], [width1;(s1+s2)*ones(size(x));width1(1)],'-', ...
     [x;flipud(x);0], [width2a;flipud(width2b);width2a(1)], '-', ...
     [x;flipud(x);0], [width3;zeros(size(x));width3(1)], '-');




%%%%%%%%%%  step responses for first solution  %%%%%%%%%%

% conductance and capacitance
G = reshape(GG*[1;sizes(:,1)],N,N);
C = reshape(CC*[1;sizes(:,1)],N,N);

% state space model
A = -inv(C)*G;
B = inv(C)*G*[ones(n,1), zeros(n,2); zeros(n,1), ones(n,1), zeros(n,1);
              zeros(n,2), ones(n,1)];
CCC = eye(N);
D = zeros(N,3);


% calculate response to step at 1st input
figure(3); hold off;

T = linspace(0,2*delays(1),1000);
[Y1,X1] = step(A,B,CCC,D,1,T);
hold off; plot(T,Y1(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y1(1000,n),'v1');
text(T(1000),Y1(1000,2*n),'v2');
text(T(1000),Y1(1000,3*n),'v3');
axis([0 2*delays(1) -0.1 1.1]); 

% show dominant time constant
plot(delays(1)*[1;1], [-0.1;1.1], '--');


% response to step at 2nd input
figure(4);  hold off;

[Y2,X2] = step(A,B,CCC,D,2,T);
hold off; plot(T,Y2(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y2(1000,n),'v1');
text(T(1000),Y2(1000,2*n),'v2');
text(T(1000),Y2(1000,3*n),'v3');
axis([0 2*delays(1) -0.1 1.1]); 

plot(delays(1)*[1;1], [-0.1;1.1], '--');
text(delays(1),-0.1,'T');


% response to step at 3rd input
figure(5);  hold off;

[Y3,X3] = step(A,B,CCC,D,3,T);
hold off; plot(T,Y3(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y3(1000,n),'v1');
text(T(1000),Y3(1000,2*n),'v2');
text(T(1000),Y3(1000,3*n),'v3');
axis([0 2*delays(1) -0.1 1.1]); 

plot(delays(1)*[1;1], [-0.1;1.1], '--');


% response to step at 1st and 2nd input
figure(6); hold off;

B = inv(C)*G*[ones(2*n,1); zeros(n,1)];
D = zeros(N,1);

T = linspace(0,2*delays(1),1000);
[Y4,X4] = step(A,B,CCC,D,1,T);
hold off; plot(T,Y4(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y4(1000,n),'v1');
text(T(1000),Y4(1000,2*n),'v2');
text(T(1000),Y4(1000,3*n),'v3');
axis([0 2*delays(1) -0.1 1.1]); 
plot(delays(1)*[1;1], [-0.1;1.1], '--');


% response to step at 1st and 3rd input
figure(7);  hold off;

B = inv(C)*G*[ones(n,1); zeros(n,1); ones(n,1)];

[Y5,X5] = step(A,B,CCC,D,1,T);
hold off; plot(T,Y5(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y5(1000,n),'v1');
text(T(1000),Y5(1000,2*n),'v2');
text(T(1000),Y5(1000,3*n),'v3');
axis([0 2*delays(1) -0.1 1.1]); 
plot(delays(1)*[1;1], [-0.1;1.1], '--');


% response to step at 2nd and 3rd wire 
figure(8);  hold off;
B = inv(C)*G*[zeros(n,1); ones(2*n,1)];

[Y6,X6] = step(A,B,CCC,D,1,T);
hold off; plot(T,Y6(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y6(1000,n),'v1');
text(T(1000),Y6(1000,2*n),'v2');
text(T(1000),Y6(1000,3*n),'v3');
axis([0 2*delays(1) -0.1 1.1]); 
plot(delays(1)*[1;1], [-0.1;1.1], '--');




%%%%%%%%%%  show second solution  %%%%%%%%%%

pause;

hold off;
colormap(gray);

w1 = sizes(1:m,2);
w2 = sizes(m+[1:m],2);
w3 = sizes(2*m+[1:m],2);
s1 = sizes(5*m+1,2);
s2 = sizes(5*m+2,2);

x = zeros(10,1);
x([1:2:2*m-1],:) = [0:m-1]';
x([2:2:2*m],:) = [1:m]';

width1 = zeros(2*m,1);
width1([1:2:2*m-1]) = s1+s2-w1;
width1([2:2:2*m]) = s1+s2-w1;

width2a = zeros(2*m,1);
width2a([1:2:2*m-1]) = s2+w2/2;
width2a([2:2:2*m]) = s2+w2/2;
width2b = zeros(2*m,1);
width2b([1:2:2*m-1]) = s2-w2/2;
width2b([2:2:2*m]) = s2-w2/2;

width3 = zeros(2*m,1);
width3([1:2:9]) = w3;
width3([2:2:10]) = w3;

plot([x;flipud(x);0], [width1;(s1+s2)*ones(size(x));width1(1)],'-', ...
     [x;flipud(x);0], [width2a;flipud(width2b);width2a(1)], '-', ...
     [x;flipud(x);0], [width3;zeros(size(x));width3(1)], '-');
axis([-0.1, m+0.1,-0.1, s1+s2+0.1]);
hold on
fill([x;m;0]',[width1;s1+s2;s1+s2]', 0.9*ones(size([x;m;0]')));
fill([x;flipud(x)]',[width2a;flipud(width2b)]', ...
     0.9*ones(size([x;x]')));
fill([x;m;0]',[width3;0;0]', 0.9*ones(size([x;0;0]')));
caxis([-1,1])
plot([x;flipud(x);0], [width1;(s1+s2)*ones(size(x));width1(1)],'-', ...
     [x;flipud(x);0], [width2a;flipud(width2b);width2a(1)], '-', ...
     [x;flipud(x);0], [width3;zeros(size(x));width3(1)], '-');




%%%%%%%%%%  step responses for second solution  %%%%%%%%%%

% conductance and capacitance
G = reshape(GG*[1;sizes(:,2)],N,N);
C = reshape(CC*[1;sizes(:,2)],N,N);

% state space model
A = -inv(C)*G;
B = inv(C)*G*[ones(n,1), zeros(n,2); zeros(n,1), ones(n,1), zeros(n,1);
              zeros(n,2), ones(n,1)];
CCC = eye(N);
D = zeros(N,3);


% calculate response to step at 1st input
figure(3); hold off;

T = linspace(0,2*delays(2),1000);
[Y1,X1] = step(A,B,CCC,D,1,T);
hold off; plot(T,Y1(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y1(1000,n),'v1');
text(T(1000),Y1(1000,2*n),'v2');
text(T(1000),Y1(1000,3*n),'v3');
axis([0 2*delays(2) -0.1 1.1]); 

% show dominant time constant
plot(delays(2)*[1;1], [-0.1;1.1], '--');


% response to step at 2nd input
figure(4);  hold off;

[Y2,X2] = step(A,B,CCC,D,2,T);
hold off; plot(T,Y2(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y2(1000,n),'v1');
text(T(1000),Y2(1000,2*n),'v2');
text(T(1000),Y2(1000,3*n),'v3');
axis([0 2*delays(2) -0.1 1.1]); 

plot(delays(2)*[1;1], [-0.1;1.1], '--');
text(delays(2),-0.1,'T');


% response to step at 3rd input
figure(5);  hold off;

[Y3,X3] = step(A,B,CCC,D,3,T);
hold off; plot(T,Y3(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y3(1000,n),'v1');
text(T(1000),Y3(1000,2*n),'v2');
text(T(1000),Y3(1000,3*n),'v3');
axis([0 2*delays(2) -0.1 1.1]); 

plot(delays(2)*[1;1], [-0.1;1.1], '--');


% response to step at 1st and 2nd input
figure(6); hold off;

B = inv(C)*G*[ones(2*n,1); zeros(n,1)];
D = zeros(N,1);

T = linspace(0,2*delays(2),1000);
[Y4,X4] = step(A,B,CCC,D,1,T);
hold off; plot(T,Y4(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y4(1000,n),'v1');
text(T(1000),Y4(1000,2*n),'v2');
text(T(1000),Y4(1000,3*n),'v3');
axis([0 2*delays(2) -0.1 1.1]); 
plot(delays(2)*[1;1], [-0.1;1.1], '--');


% response to step at 1st and 3rd input
figure(7);  hold off;

B = inv(C)*G*[ones(n,1); zeros(n,1); ones(n,1)];

[Y5,X5] = step(A,B,CCC,D,1,T);
hold off; plot(T,Y5(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y5(1000,n),'v1');
text(T(1000),Y5(1000,2*n),'v2');
text(T(1000),Y5(1000,3*n),'v3');
axis([0 2*delays(2) -0.1 1.1]); 
plot(delays(2)*[1;1], [-0.1;1.1], '--');


% response to step at 2nd and 3rd wire 
figure(8);  hold off;
B = inv(C)*G*[zeros(n,1); ones(2*n,1)];

[Y6,X6] = step(A,B,CCC,D,1,T);
hold off; plot(T,Y6(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y6(1000,n),'v1');
text(T(1000),Y6(1000,2*n),'v2');
text(T(1000),Y6(1000,3*n),'v3');
axis([0 2*delays(2) -0.1 1.1]); 
plot(delays(2)*[1;1], [-0.1;1.1], '--');

