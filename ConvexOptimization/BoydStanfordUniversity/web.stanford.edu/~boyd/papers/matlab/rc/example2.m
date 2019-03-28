% 
% EXAMPLE 2:  Wire and driver sizing
% 


% circuit parameters

n = 21;        % number of nodes per wire 
m = n-1;       % number of segments per wire
g = 1.0;       % output conductance is g times driver size 
c0 = 1.0;      % input capacitance of driver is co + c*driver size 
c = 3.0; 
alpha = 10;    % wire segment: two capacitances beta*width
beta = 0.5;    % wire segment: conductance alpha*width
C = 50;        % external load
L = 10.0;      % area is sum xi + L*(d1+d2)
wmax = 2.0;    % maximum wire width
dmax = 100.0;  % maximum driver size 


% 2*m+2 variables x
% x(1:m):      widths of first wire
% x(m+1:2*m):  widths of 2nd wire
% x(2*m+1):    size of 1st driver 
% x(2*m+2):    size of 2nd driver


% capacitance matrix of first circuit: 
CC1 = zeros(n*n,2*m+3);

%  external load from second driver
CC1(n*n,1) = c0;
CC1(n*n,2*m+3) = c;

% capacitances from segments
for i=1:(n-1);
   CC1((i-1)*n+i,i+1) = CC1((i-1)*n+i,i+1) + beta;
   CC1(i*n+i+1,i+1) = CC1(i*n+i+1,i+1) + beta;
end;

% conductance matrix of first circuit
GG1 = zeros(n*n,2*m+3);

% (1,1) element of GG1_0 is g*d1
GG1(1,2*m+2) = g;

% conductances from segments
for i=1:(n-1)         
   GG1((i-1)*n+i,i+1) = GG1((i-1)*n+i,i+1) + alpha;
   GG1((i-1)*n+i+1,i+1) = GG1((i-1)*n+i+1,i+1) - alpha;
   GG1(i*n+i,i+1) = GG1(i*n+i,i+1) - alpha;
   GG1(i*n+i+1,i+1) = GG1(i*n+i+1,i+1) + alpha;
end;

% capacitance matrix of second circuit 
CC2 = zeros(n*n,2*m+3);

% external load
CC2(n*n,1) = C;

% capacitances from segments
for i=1:(n-1);
   CC2((i-1)*n+i,m+i+1) = CC2((i-1)*n+i,m+i+1)  + beta;
   CC2(i*n+i+1,m+i+1) = CC2(i*n+i+1,m+i+1)  + beta;
end;


% conductance matrix of second circuit
GG2 = zeros(n*n, 2*m+3);

% conductance of second driver
GG2(1,2*m+3) = g;

% conductances of segments
for i=1:(n-1);         
   GG2((i-1)*n+i,m+i+1) = GG2((i-1)*n+i,m+i+1) + alpha;
   GG2((i-1)*n+i+1,m+i+1) = GG2((i-1)*n+i+1,m+i+1) - alpha;
   GG2(i*n+i,m+i+1) = GG2(i*n+i,m+i+1) - alpha;
   GG2(i*n+i+1,m+i+1) = GG2(i*n+i+1,m+i+1) + alpha;
end;


% LMI-constraints:
% (-1/delay)*C1 + G1 (+t*I) >= 0   (nxn)
% (-1/delay)*C2 + G2 (+t*I) >= 0   (nxn)
% x1 >= 0                          (m scalar constr.)
% x2 >= 0                          (m scalar constr.)
% d1 >= 0                          (scalar)
% d2 >= 0                          (scalar)
% x1 <= wmax                       (m scalar constr.)
% x2 <= wmax                       (m scalar constr.)
% d1 <= dmax                       (scalar)
% d2 <= dmax                       (scalar)
%
% variables:
% x1 (m) 
% x2 (m)
% d1 (scalar)
% d2 (scalar)
% t  (scalar), only for phase 1

% store LMI in matrix FF
FF = zeros(2*n*n+4*m+4, 2*m+4); 
blck_szs = [n; n; ones(4*m+4,1)];

% fill in first two blocks later

% last column is identity (used only during phase 1)
FF(1:n*n,2*m+4) = reshape(eye(n),n*n,1);
FF(n*n+[1:n*n],2*m+4) = reshape(eye(n),n*n,1);

% lower bounds 
FF(2*n*n+[1:2*m+2],1+[1:2*m+2]) = eye(2*m+2);

% upper bounds 
FF(2*n*n+2*m+2+[1:2*m+2],1+[1:2*m+2]) = -eye(2*m+2);
FF(2*n*n+2*m+2+[1:2*m+2],1) = [wmax*ones(2*m,1); dmax; dmax];

% cost function
c = [ones(2*m,1); L; L ];  


%%%%%%%%%%  compute tradeoff curve  %%%%%%%%%%

% tolerances for SP
abstol = 1e-4;
reltol = 1e-4;

% points on the tradoff curve
nopts = 50;
delays = linspace(150,500,nopts);
areas = zeros(1,nopts);

% compute tradeoff curve
for i=1:nopts

   delay = delays(i);

   disp([' ']);
   disp(['Point ', int2str(i), ...
        ' on the tradeoff curve (delay = ', num2str(delay), ').']);

   % first block of F(x) = G1(x) - (2/delay)*C1(x)
   FF(1:n*n,1:2*m+3) = GG1 - (2/delay)*CC1;

   % second block of F(x) = G2(x) - (2/delay)*C2(x)
   FF(n*n+[1:n*n],1:2*m+3) = GG2 - (2/delay)*CC2;

   % try very thick wire, very large drivers
   x0 = 0.99*[wmax*ones(2*m,1); dmax; dmax];
   t = min(min(eig(reshape(FF(1:n*n,1:2*m+3)*[1;x0],n,n))), ...
           min(eig(reshape(FF(2*n+[1:n*n],1:2*m+3)*[1;x0],n,n))));

   feasible = 1;
   if (t <= 0)   % x0 is infeasible, run phase 1
   
      disp([' ']);  disp(['Phase 1.']);

      % use slightly tighter bounds on x1, x2, d1, d2 to ensure 
      % strict feasibility
      FF(2*n*n+[1:2*m+2],1) = -1e-5*ones(2*m+2,1);
      FF(2*n*n+2*m+2+[1:2*m+2],1) = ...
         (1-1e-5)*[ wmax*ones(2*m,1); dmax; dmax ];
    
      % dual feasible solution
      z = FF(1:2*n*n,2:2*m+3)' * ...
          [ reshape(eye(n)/(2*n),n*n,1); reshape(eye(n)/(2*n),n*n,1) ];
      Z0 = [ reshape(eye(n)/(2*n),n*n,1); reshape(eye(n)/(2*n),n*n,1); 
             -min(z,0)+1e-3; max(z,0)+1e-3 ]; 

      % call SP
      [x,Z,ul,info,time] = sp(FF, blck_szs, [zeros(2*m+2,1);1], ...
          [x0; -1.1*t], Z0, 10.0, abstol, -1.0, 0.0, 100);

      if (ul(1) >= 0),
         areas(i) = Inf;
         if (ul(2) >= 0),
            disp(['Infeasible.']);
         else
            disp(['Feasibility could not be determined.']);
         end;
         feasible = 0;
      else x0 = x(1:2*m+2);    % initial point for phase 2.
      end;

   end;

   if (feasible) % use x0 as starting point, run phase 2

      disp([' ']); disp(['Phase 2.']);

      % use original bounds on x again
      FF(2*n*n+[1:2*m+2],1) = zeros(2*m+2,1);
      FF(2*n*n+2*m+2+[1:2*m+2],1) = [wmax*ones(2*m,1); dmax; dmax ];

      % dual feasible point
      z = FF(1:2*n*n,2:2*m+3)' * ...
         [ reshape(eye(n)/(2*n),n*n,1); ...
           reshape(eye(n)/(2*n),n*n,1) ] - c;
      Z0 = [ reshape(eye(n)/(2*n),n*n,1); 
             reshape(eye(n)/(2*n),n*n,1); 
             -min(z,0)+1e-3; max(z,0)+1e-3 ]; 

      % call SP
      [x,Z,ul,info,time] = sp(FF(:,1:2*m+3), blck_szs, c, ...
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



%%%%%%%%%%  draw one solution  %%%%%%%%%% 

delay = 189;

% tolerances for SP
abstol = 1e-6;
reltol = 1e-6;

disp([' ']);
disp(['Compute solution for delay = ', num2str(delay), '.']);

% first block of F(x) = G1(x) - (2/delay)*C1(x)
FF(1:n*n,1:2*m+3) = GG1 - (2/delay)*CC1;

% second block of F(x) = G2(x) - (2/delay)*C2(x)
FF(n*n+[1:n*n],1:2*m+3) = GG2 - (2/delay)*CC2;

% try very thick wire, very large drivers
x0 = 0.99*[wmax*ones(2*m,1); dmax; dmax];
t = min(min(eig(reshape(FF(1:n*n,1:2*m+3)*[1;x0],n,n))), ...
        min(eig(reshape(FF(2*n+[1:n*n],1:2*m+3)*[1;x0],n,n))));

feasible = 1;
if (t <= 0)   % x0 is infeasible, run phase 1
   
   disp([' ']);  disp(['Phase 1.']);

   % use slightly tighter bounds on x1, x2, d1, d2 to ensure 
   % strict feasibility
   FF(2*n*n+[1:2*m+2],1) = -1e-5*ones(2*m+2,1);
   FF(2*n*n+2*m+2+[1:2*m+2],1) = ...
      (1-1e-5)*[ wmax*ones(2*m,1); dmax; dmax ];
    
   % dual feasible solution
   z = FF(1:2*n*n,2:2*m+3)' * ...
       [ reshape(eye(n)/(2*n),n*n,1); reshape(eye(n)/(2*n),n*n,1) ];
   Z0 = [ reshape(eye(n)/(2*n),n*n,1); reshape(eye(n)/(2*n),n*n,1); 
          -min(z,0)+1e-3; max(z,0)+1e-3 ]; 

   % call SP
   [x,Z,ul,info,time] = sp(FF, blck_szs, [zeros(2*m+2,1);1], ...
       [x0; -1.1*t], Z0, 10.0, abstol, -1.0, 0.0, 100);

   if (ul(1) >= 0),
      areas(i) = Inf;
      if (ul(2) >= 0),
         disp(['Infeasible.']);
      else
         disp(['Feasibility could not be determined.']);
      end;
      feasible = 0;
   else x0 = x(1:2*m+2);    % initial point for phase 2.
   end;

end;

if (feasible) % use x0 as starting point, run phase 2

   disp([' ']); disp(['Phase 2.']);

   % use original bounds on x again
   FF(2*n*n+[1:2*m+2],1) = zeros(2*m+2,1);
   FF(2*n*n+2*m+2+[1:2*m+2],1) = [wmax*ones(2*m,1); dmax; dmax ];

   % dual feasible point
   z = FF(1:2*n*n,2:2*m+3)' * ...
      [ reshape(eye(n)/(2*n),n*n,1); ...
        reshape(eye(n)/(2*n),n*n,1) ] - c;
   Z0 = [ reshape(eye(n)/(2*n),n*n,1); 
          reshape(eye(n)/(2*n),n*n,1); 
          -min(z,0)+1e-3; max(z,0)+1e-3 ]; 

   % call SP
   [x,Z,ul,info,time] = sp(FF(:,1:2*m+3), blck_szs, c, ...
       x0, Z0, 10.0, abstol, reltol, 0.0, 100);

end;


figure(2)   
colormap(gray)

os = 3;  % plot drivers as block with width os and height L/(2*os);

% contour of  first wire
width1 = zeros(2*m,1);
width1([1:2:2*m-1],:) = x(1:m);
width1([2:2:2*m],:) = x(1:m);
x1 = zeros(2*m,1);
x1([1:2:2*m-1],:) = [0:m-1]';
x1([2:2:2*m],:) = [1:m]';

% contour of second wire
width2 = zeros(2*m,1);
width2([1:2:2*m-1],:) = x(m+[1:m]);
width2([2:2:2*m],:) = x(m+[1:m]);
x2 = zeros(2*m,1);
x2([1:2:2*m-1],:) = [m:2*m-1]';
x2([2:2:2*m],:) = [m+1:2*m]';

% drivers
d1=x(2*m+1);  
d2=x(2*m+2);
xd1=[0;os;os;0;0];  
xd2=[os+m+0.5; os+m+0.5+os; os+m+0.5+os; os+m+0.5; os+m+0.5];


% draw four contours with 0.5 space between them
hold off
plot((os+0.5)+[x1;flipud(x1);0], ...
     [-0.5*width1; flipud(0.5*width1); -0.5*width1(1)], '-',  ...
     (os+1.5+os)+[x2;flipud(x2);m], ...
     [-0.5*width2; flipud(0.5*width2); -0.5*width2(1)],'-', ...
     xd1, (L/(2*os))*[-d1;-d1;d1;d1;-d1], '-', ...
     0.5+xd2, (L/(2*os))*[-d2;-d2;d2;d2;-d2], '-'); 
hold on;

% shade
fill((os+0.5)+[x1;flipud(x1);0], ...
     [-0.5*width1; flipud(0.5*width1); -0.5*width1(1)], ... 
     0.9*ones(size([x1;x1;0]')));
fill((os+1.5+os)+[x2;flipud(x2);m], ...
     [-0.5*width2; flipud(0.5*width2); -0.5*width2(1)], ...
     0.9*ones(size([x2;x2;0]')));
fill(xd1, (L/(2*os))*[-d1;-d1;d1;d1;-d1], 0.8*ones(size(xd1)));
fill(xd2+0.5,(L/(2*os))*[-d2;-d2;d2;d2;-d2], 0.8*ones(size(xd2))); 
plot((os+0.5)+[x1;flipud(x1);0], ...
         [-0.5*width1; flipud(0.5*width1); -0.5*width1(1)], '-',  ...
     (os+1.5+os)+[x2;flipud(x2);m], ...
         [-0.5*width2; flipud(0.5*width2); -0.5*width2(1)],'-', ...
     xd1, (L/(2*os))*[-d1;-d1;d1;d1;-d1], '-', ...
     xd2+0.5, (L/(2*os))*[-d2;-d2;d2;d2;-d2], '-'); 
caxis([-1,1]);

% reset axes
axis([-0.1 os+os+2*m+1.6 -5 5]);
set(gca,'XTick',[os+0.5 os+ceil(m/2)+.5 os+m+0.5 ...
    2*os+m+1.5 2*os+m+ceil(m/2)+1.5 2*os+2*m+1.5]);
set(gca,'XTickLabels',[' 0'; sprintf('%2d',ceil(m/2)); ...
    sprintf('%2d',m); ' 0'; sprintf('%2d',ceil(m/2)); ...
    sprintf('%2d',m)]);




%%%%%%%%%%  step responses  %%%%%%%%%%

figure(3); 
hold off;

% first circuit
C1 = reshape(CC1*[1;x],n,n);
G1 = reshape(GG1*[1;x],n,n);

% linear system:  C1*vdot = -G1*v + G1*e*u;   y = v 
A1 = -inv(C1)*G1;   
B1 = inv(C1)*G1*ones(n,1);
CCC1 = eye(n);
D1 = zeros(n,1);

% dominant time constant and Elmore delay first wire
tdom1=max(eig(inv(G1)*C1));
telm1=max(sum((inv(G1)*C1)'));

% second circuit
C2 = reshape(CC2*[1;x],n,n);
G2 = reshape(GG2*[1;x],n,n);

% linear system:  C2*vdot = -G2*v + G2*e*u;   y = v 
A2 = -inv(C2)*G2;   
B2 = inv(C2)*G2*ones(n,1);
CCC2 = eye(n);
D2 = zeros(n,1);

% dominant time constant and Elmore delay second wire
tdom2=max(eig(inv(G2)*C2));
telm2=max(sum((inv(G2)*C2)'));

% 2 step responses
T = linspace(0,1000,1000);
[Y1,X1] = step(A1,B1,CCC1,D1,1,T);
[Y2,X2] = step(A2,B2,CCC2,D2,1,T);

% 50% threshold delays
ind = find(Y1(:,n) >= 0.5);  tthres1 = ind(1);
ind = find(Y2(:,n) >= 0.5);  tthres2 = ind(1);

% plot step responses
plot(T,Y1(:,n),'-', T+tdom1,Y2(:,n),'-');  hold on;
axis([0 T(500) 0 1]);
xlabel('time');
ylabel('v');

% various delays
text(tdom1,0,'d1');
text(telm1,0,'e1');
text(tthres1,0,'t1');
text(tdom1+tdom2,0,'d2');
text(tdom1+telm2,0,'e2');
text(tdom1+tthres2,0,'t2');
plot(tdom1*[1;1],[0;1],'--');
plot(telm1*[1;1],[0;1],'--');
plot(tthres1*[1;1],[0;1],'--');
plot((tdom1+tdom2)*[1;1],[0;1],'--');
plot((tdom1+telm2)*[1;1],[0;1],'--');
plot((tdom1+tthres2)*[1;1],[0;1],'--');

