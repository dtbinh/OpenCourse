%
% EXAMPLE 3:  Wire sizing and topology design
%

% circuit parameters
n = 4;             % number of nodes
m = 6;             % number of branches
G = 0.1;           % resistor between node 1 and 0
Co = 10;           % load capacitance
alpha1 = 1.0;      % conductance per segment
beta1 = 10;        % capacitance per segment is 2*beta
alpha2 = 1.0;      
beta2 = 10;
alpha3 = 1.0;     
beta3 = 100;
alpha4 = 1.0;    
beta4 = 1;
alpha5 = 1.0;   
beta5 = 1;
alpha6 = 1.0;  
beta6 = 1;
wmax = 10.0;  


% capacitance and conductance matrix 
CC = zeros(n*n,m+1);
GG = zeros(n*n,m+1);

% constant terms
CC(2*n+3,1) = Co;
GG(1,1) = G;

% branch 1;
CC(1,2) = CC(1,2) + beta1;        
CC(n+2,2) = CC(n+2,2) + beta1;    
GG(1,2) = GG(1,2) + alpha1;
GG(2,2) = GG(2,2) - alpha1; 
GG(n+1,2) = GG(n+1,2) - alpha1;
GG(n+2,2) = GG(n+2,2) + alpha1;

% branch 2;
CC(n+2,3) = CC(n+2,3) + beta2;        
CC(2*n+3,3) = CC(2*n+3,3) + beta2;   
GG(n+2,3) = GG(n+2,3) + alpha2;
GG(n+3,3) = GG(n+3,3) - alpha2; 
GG(2*n+2,3) = GG(2*n+2,3) - alpha2;
GG(2*n+3,3) = GG(2*n+3,3) + alpha2;

% branch 3;
CC(1,4) = CC(1,4) + beta3;         
CC(2*n+3,4) = CC(2*n+3,4) + beta3;    
GG(1,4) = GG(1,4) + alpha3;
GG(3,4) = GG(3,4) - alpha3; 
GG(2*n+1,4) = GG(2*n+1,4) - alpha3;
GG(2*n+3,4) = GG(2*n+3,4) + alpha3;

% branch 4;
CC(1,5) = CC(1,5) + beta4;         
CC(3*n+4,5) = CC(3*n+4,5) + beta4;    
GG(1,5) = GG(1,5) + alpha4;
GG(4,5) = GG(5,5) - alpha4; 
GG(3*n+1,5) = GG(3*n+1,5) - alpha4;
GG(3*n+4,5) = GG(3*n+4,5) + alpha4;

% branch 5;
CC(n+2,6) = CC(n+2,6) + beta5;         
CC(3*n+4,6) = CC(3*n+4,6) + beta5;    
GG(n+2,6) = GG(n+2,6) + alpha5;
GG(n+4,6) = GG(n+4,6) - alpha5; 
GG(3*n+2,6) = GG(3*n+2,6) - alpha5;
GG(3*n+4,6) = GG(3*n+4,6) + alpha5;

% branch 6;
CC(3*n+4,7) = CC(3*n+4,7) + beta6;         
CC(2*n+3,7) = CC(2*n+3,7) + beta6;    
GG(3*n+4,7) = GG(3*n+4,7) + alpha6;
GG(3*n+3,7) = GG(3*n+3,7) - alpha6; 
GG(2*n+4,7) = GG(2*n+4,7) - alpha6;
GG(2*n+3,7) = GG(2*n+3,7) + alpha6;

% 
% LMI constraints:
% (-1/delay)*C + G  (+t*I) >= 0  (nxn)
% x >= 0                         (m scalar constraints)
% x <= wmax                      (m scalar constraints)

FF = zeros(n*n+2*m, m+2);    % last column used in phase 1 only
blck_szs = [n; ones(2*m,1)];

% fill in 1st block later

% last column is identity
FF(1:n*n,m+2) = reshape(eye(n),n*n,1);

% lower and upper bounds
FF(n*n+[1:m],2:m+1) = eye(m);
FF(n*n+m+[1:m],2:m+1) = -eye(m);
FF(n*n+m+[1:m],1) = wmax*ones(m,1);

% cost function
c = ones(m,1);


%%%%%%%%%%  compute tradeoff curve  %%%%%%%%%%

% tolerances for SP
abstol = 1e-4;
reltol = 1e-4;

% points on tradeoff curve
nopts = 50;
delays = linspace(180,800,nopts);
areas = zeros(1,nopts);

% compute tradeoff curve
for i=1:nopts
   
   delay = delays(i);
   disp([' ']);
   disp(['Point ', int2str(i), ...
        ' on the tradeoff curve (delay = ', num2str(delay), ').']);

   % first block of F(x) = G(x) - (1/delay)*C(x)
   FF(1:n*n,1:m+1) = GG - (1/delay)*CC;

   % try very thick wire
   x0 = 0.99*wmax*ones(m,1);
   t = min(min(eig(reshape(FF(1:n*n,1:m+1)*[1;x0],n,n))));

   feasible = 1;
   if (t <= 0)  % x0 is infeasible, run phase 1

      disp([' ']);  disp(['Phase 1.']);

      % use slightly tighter bounds on x to ensure strict feasibility
      FF(n*n+[1:m],1) = -1e-5*ones(m,1);
      FF(n*n+m+[1:m],1) = (1-1e-5)*wmax*ones(m,1);

      % dual feasible solution
      z = FF(1:n*n,2:m+1)' * reshape(eye(n)/n,n*n,1);
      Z0 = [ reshape(eye(n)/n,n*n,1); -min(z,0)+1e-3; max(z,0)+1e-3 ];

      % call SP
      [x,Z,ul,info,time] = sp(FF, blck_szs, [zeros(m,1);1], ...
          [x0; -1.1*t], Z0, 10.0, abstol, -1.0, 0.0, 100);

      if (ul(1) >= 0),
         areas(i) = Inf;
         if (ul(2) >= 0),
            disp(['Infeasible.']);
         else
            disp(['Feasibility could not be determined.']);
         end;
         feasible = 0;
      else x0 = x(1:m);    % initial point for phase 2.
      end;

   end;

   if (feasible) % use x0 as starting point, run phase 2

      disp([' ']); disp(['Phase 2.']);

      % use original bounds on x again
      FF(n*n+[1:m],1) = zeros(m,1);
      FF(n*n+m+[1:m],1) = wmax*ones(m,1); 

      % dual feasible point
      z = FF(1:n*n,2:m+1)' * reshape(eye(n)/n,n*n,1) - c;
      Z0 = [ reshape(eye(n)/n,n*n,1); -min(z,0)+1e-3; max(z,0)+1e-3 ];

      % call SP
      [x,Z,ul,info,time] = sp(FF(:,1:m+1), blck_szs, c, ...
          x0, Z0, 10.0, abstol, reltol, 0.0, 100);

      % optimal value
      areas(i) = c'*x;

   end;

end;

figure(1);
ind = finite(areas);
plot(areas(ind), delays(ind));
xlabel('area');
ylabel('Tdom');


%%%%%%%%%%  compute three solutions  %%%%%%%%%%


% tolerances for SP
abstol = 1e-8;
reltol = 1e-8;

delays = [200 400 600];
sizes = zeros(m,3);

for i=[1,2,3]
   delay = delays(i);
   disp([' ']);
   disp(['Compute solution ', int2str(i), ...
        ' (delay = ', num2str(delay), ').']);

   % first block of F(x) = G(x) - (1/delay)*C(x)
   FF(1:n*n,1:m+1) = GG - (1/delay)*CC;

   % try very thick wire
   x0 = 0.99*wmax*ones(m,1);
   t = min(min(eig(reshape(FF(1:n*n,1:m+1)*[1;x0],n,n))));

   feasible = 1;
   if (t <= 0)  % x0 is infeasible, run phase 1

      disp([' ']);  disp(['Phase 1.']);

      % use slightly tighter bounds on x to ensure strict feasibility
      FF(n*n+[1:m],1) = -1e-5*ones(m,1);
      FF(n*n+m+[1:m],1) = (1-1e-5)*wmax*ones(m,1);

      % dual feasible solution
      z = FF(1:n*n,2:m+1)' * reshape(eye(n)/n,n*n,1);
      Z0 = [ reshape(eye(n)/n,n*n,1); -min(z,0)+1e-3; max(z,0)+1e-3 ];

      % call SP
      [x,Z,ul,info,time] = sp(FF, blck_szs, [zeros(m,1);1], ...
          [x0; -1.1*t], Z0, 10.0, abstol, -1.0, 0.0, 100);

      if (ul(1) >= 0),
         areas(i) = Inf;
         if (ul(2) >= 0),
            disp(['Infeasible.']);
         else
            disp(['Feasibility could not be determined.']);
         end;
         feasible = 0;
      else x0 = x(1:m);    % initial point for phase 2.
      end;

   end;

   if (feasible) % use x0 as starting point, run phase 2

      disp([' ']); disp(['Phase 2.']);

      % use original bounds on x again
      FF(n*n+[1:m],1) = zeros(m,1);
      FF(n*n+m+[1:m],1) = wmax*ones(m,1); 

      % dual feasible point
      z = FF(1:n*n,2:m+1)' * reshape(eye(n)/n,n*n,1) - c;
      Z0 = [ reshape(eye(n)/n,n*n,1); -min(z,0)+1e-3; max(z,0)+1e-3 ];

      % call SP
      [sizes(:,i),Z,ul,info,time] = sp(FF(:,1:m+1), blck_szs, c, ...
          x0, Z0, 10.0, abstol, reltol, 0.0, 100);

   end;


end;

disp(['Three solutions:']);
sizes


%%%%%%%%%%  step responses  %%%%%%%%%%

figure(3)

hold off
x = sizes(:,1);

% conductance, capacitance matrix
G = reshape(GG*[1;x],n,n);
C = reshape(CC*[1;x],n,n);

% linear system:  C*vdot = -G*v + G*e*u;   y = v
A = -inv(C)*G;
B = inv(C)*G*ones(n,1);
CCC = eye(n);
D = zeros(4,1);

% step response
T = linspace(0,1000,1000);
[Y,X] = step(A,B,CCC,D,1,T);
hold off; plot(T,Y(:,[1,3,4]),'-');  hold on;

% calculate 3 delays
ind = find(Y(:,3) >= 0.5);
thresh = T(ind(1));
tdom=max(eig(inv(G)*C));
elmore=max(sum((inv(G)*C)'));
plot(tdom*[1;1], [0;1],'--', elmore*[1;1], [0;1],'--', ...
     thresh*[1;1], [0;1], '--');
text(tdom,0,'d');
text(elmore,0,'e');
text(thresh,0,'t');
title('step responese for solution a');


figure(4)
hold off
x = sizes(:,2);

% conductance, capacitance matrix
G = reshape(GG*[1;x],n,n);
C = reshape(CC*[1;x],n,n);

% linear system:  C*vdot = -G*v + G*e*u;   y = v
A = -inv(C)*G;
B = inv(C)*G*ones(n,1);
CCC = eye(n);
D = zeros(4,1);

% step response
T = linspace(0,1000,1000);
[Y,X] = step(A,B,CCC,D,1,T);
hold off; plot(T,Y(:,[1,3,4]),'-');  hold on;

% calculate 3 delays
ind = find(Y(:,3) >= 0.5);
thresh = T(ind(1));
tdom=max(eig(inv(G)*C));
elmore=max(sum((inv(G)*C)'));
plot(tdom*[1;1], [0;1],'--', elmore*[1;1], [0;1],'--', ...
     thresh*[1;1], [0;1], '--');
text(tdom,0,'d');
text(elmore,0,'e');
text(thresh,0,'t');
title('step response for solution b');



figure(5)
hold off
x = sizes(:,3);

% conductance, capacitance matrix
G = reshape(GG*[1;x],n,n);
C = reshape(CC*[1;x],n,n);

% linear system:  C*vdot = -G*v + G*e*u;   y = v
A = -inv(C)*G;
B = inv(C)*G*ones(n,1);
CCC = eye(n);
D = zeros(4,1);

% step response
T = linspace(0,1000,1000);
[Y,X] = step(A,B,CCC,D,1,T);
hold off; plot(T,Y(:,[1,3]),'-');  hold on;

% calculate 3 delays
ind = find(Y(:,3) >= 0.5);
thresh = T(ind(1));
tdom=max(eig(inv(G)*C));
elmore=max(sum((inv(G)*C)'));
plot(tdom*[1;1], [0;1],'--', elmore*[1;1], [0;1],'--', ...
     thresh*[1;1], [0;1], '--');
text(tdom,0,'d');
text(elmore,0,'e');
text(thresh,0,'t');
title('step response for solution c');

