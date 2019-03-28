%
% EXAMPLE 1:  Wire sizing
%


% circuit parameters
n=21;          % number of nodes; n-1 is number of segments in the wire
m=n-1;         % number of variables = number of segments
beta = 0.5;    % segment has two capacitances beta*xi
alpha = 1;     % conductance is alpha*xi per segment
Go = 1;        % driver conductance
Co = 10;       % load capacitance
wmax = 1.0;    % upper bound on x 


%
% capacitance matrix: C(x)  = C0 + x1*C1 + x2*C2 + ... + xn*Cn 
% store Ci's in matrix CC = [ C0(:) C1(:) C2(:) ... Cn(:) ]  
% ie, CC is such that C(x) =  reshape(CC*[1;x],n,n) 
% 

CC = zeros(n*n,m+1);

% (n,n)-element of C0 is Co
CC(n*n,1) = Co;

% wire capacitances:
% (i,i) element of Ci is beta
% (i+1,i+1) element of C_i is beta
for i=1:(n-1)
   CC((i-1)*n+i,i+1) = CC((i-1)*n+i,i+1) + beta; 
   CC(i*n+i+1,i+1) = CC(i*n+i+1,i+1) + beta;
end;


% 
% conductance matrix: G(x)  = G0 + x1*G1 + x2*G2 + ... + xn*Gn 
% store Gi's in matrix GG = [ G0(:) G1(:) G2(:) ... Gn(:) ]  
% ie, GG is such that G(x) = reshape(GG*[1;x],n,n)  
% 

GG = zeros(n*n,m+1);

% (1,1) element of G0 is Go
GG(1,1) = Go;

% wire conductances
% (i,i) element of Gi is alpha
% (i+1,i+1) element of Gi is -alpha
% (i,i+1) element of Gi is -alpha
% (i+1,i+1) element of Gi is alpha
for i=1:(n-1)
   GG((i-1)*n+i,i+1) = GG((i-1)*n+i,i+1) + alpha;
   GG((i-1)*n+i+1,i+1) = GG((i-1)*n+i+1,i+1) - alpha;
   GG(i*n+i,i+1) = GG(i*n+i,i+1) - alpha;
   GG(i*n+i+1,i+1) = GG(i*n+i+1,i+1) + alpha;
end;


%
%
% PHASE 1
%
% primal SDP:  minimize   t
%              subject to (-1/delay)*C(x) + G(x) + t*I >= 0
%                         1e-5 <= x <= (1-1e-5)*wmax
%
% dual SDP:    maximize   Tr ((1/delay)*C0 - G0)*Z - 1e-5*e^T*zl
%                            + (1-1e-5)*wmax*e^T*zu  
%              subject to Tr ((-1/delay)*Ci + Gi)*Z - zu + zl = 0
%                         Tr Z = 1
%                         Z >= 0, zu >= 0, zl >= 0
%
% PHASE 2:
%
% primal SDP:  minimize   sum xi
%              subject to (-1/delay)*C(x) + G(x) >= 0
%                         0 <= x <= wmax
% 
% dual SDP:    maximize   Tr ((1/delay)*C0 - G0)*Z + wmax*e^T*zu
%              subject to Tr ((-1/delay)*Ci + Gi)*Z - zu + zl = e
%                         Z >= 0, zu >= 0, zl >= 0
% 
 
% store LMI in matrix FF
FF = zeros(n*n+2*m,m+2);   % last column used only in phase 1
blck_szs = [n; ones(m,1); ones(m,1)]; 

% will fill in GG and CC later

% last column is identity (used only during phase 1)
FF(1:n*n,m+2) = reshape(eye(n),n*n,1);

% lower bounds on x;
FF(n*n+[1:m],2:m+1) = eye(m);

% upper bounds on x;
FF(n*n+m+[1:m],2:m+1) = -eye(m);
FF(n*n+m+[1:m],1) = wmax*ones(m,1);



%%%%%%%%%%  compute tradeoff curve  %%%%%%%%%% 

% tolerances for SP
abstol = 1e-4;  
reltol = 1e-4;

% points on the tradeoff curve 
nopts = 50;                          % no of points
delays = linspace(400,2000,nopts);   % begin and end point 
areas = zeros(1,nopts);

% compute tradeoff curve
for i=1:nopts

   delay = delays(i);

   disp([' ']);
   disp(['Point ', int2str(i), ...
        ' on the tradeoff curve (delay = ', num2str(delay), ').']);
 
   % first block of F(x) is G(x) - (1/delay)*C(x)
   FF(1:n*n,1:m+1) = GG - (1/delay)*CC;  

   % see if very thick wire is feasible
   x0 = 0.99*wmax*ones(m,1);
   t = min(eig(reshape(FF(1:n*n,1:m+1)*[1;x0],n,n)));

   feasible = 1;
   if (t <= 0)  % x0 is infeasible, run phase 1
 
      disp([' ']); disp(['Phase 1.']);

      % use slightly tighter bounds on x to have strict feasibility
      FF(n*n+[1:m],1) = -1e-5*ones(m,1);
      FF(n*n+m+[1:m],1) = (1-1e-5)*wmax*ones(m,1);

      % dual initial solution for phase 1
      z = FF(1:n*n,2:m+1)'*reshape(eye(n)/n,n*n,1); 
      Z0 = [ reshape(eye(n)/n,n*n,1); -min(z,0)+1e-3;  max(z,0)+1e-3];

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

      % dual initial point
      z = FF(1:n*n,2:m+1)'*reshape(eye(n)/n,n*n,1) - 1; 
      Z0 = [ reshape(eye(n)/n,n*n,1); -min(z,0)+1e-3;  max(z,0)+1e-3];

      % call SP
      [x,Z,ul,info,time] = sp(FF(:,1:m+1), blck_szs, ones(m,1), ...
          x0, Z0, 10.0, abstol, reltol, 0.0, 100);

      % optimal value 
      areas(i) = sum(x);

   end;

end;
  
figure(1)
ind = finite(areas);   % skip infeasible values
plot(areas(ind), delays(ind));
xlabel('area');
ylabel('Tdom');


%%%%%%%%%%  compute four solutions  %%%%%%%%%%

% tolerances for SP
abstol=1e-6;
reltol=1e-6;

% compute solution for these four values of delay
delays = [370, 400, 600, 1800];
sizes = zeros(m,4);

for i=1:4
 
   delay = delays(i);

   disp([' ']);
   disp(['Compute solution ', int2str(i), ...
         ' (delay = ', num2str(delay), ').']);
 
   % first block of F(x) is G(x) - (1/delay)*C(x)
   FF(1:n*n,1:m+1) = GG - (1/delay)*CC;  

   % see if very thick wire is feasible
   x0 = 0.99*wmax*ones(m,1);
   t = min(eig(reshape(FF(1:n*n,1:m+1)*[1;x0],n,n)));

   feasible = 1;
   if (t <= 0)  % x0 is infeasible, run phase 1

      disp([' ']); disp(['Phase 1.']);

      % use slightly tighter bounds on x to have strict feasibility
      FF(n*n+[1:m],1) = -1e-5*ones(m,1);
      FF(n*n+m+[1:m],1) = (1-1e-5)*wmax*ones(m,1);
 
      % dual initial solution for phase 1
      z = FF(1:n*n,2:m+1)'*reshape(eye(n)/n,n*n,1); 
      Z0 = [ reshape(eye(n)/n,n*n,1); -min(z,0)+1e-3;  max(z,0)+1e-3];

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
      else x0 = x(1:m);  % initial point for phase 2
      end;

   end;

   if (feasible)  % use x0 as starting point, run phase 2

      disp([' ']); disp(['Phase 2.']);

      % use original bounds on x again
      FF(n*n+[1:m],1) = zeros(m,1);
      FF(n*n+m+[1:m],1) = wmax*ones(m,1);

      % dual initial point
      z = FF(1:n*n,2:m+1)'*reshape(eye(n)/n,n*n,1) - 1; 
      Z0 = [ reshape(eye(n)/n,n*n,1); -min(z,0)+1e-3;  max(z,0)+1e-3];

      % call SP
      [sizes(:,i),Z,ul,info,time] = sp(FF(:,1:m+1), ...
          blck_szs, ones(m,1), x0, Z0, 10.0, abstol, reltol, ...
          0.0, 100);

   end;

end;



%%%%%%%%%%  draw wires for four solutions  %%%%%%%%%%

figure(2)   
colormap(gray);

width = zeros(2*m,4);
width([1:2:2*m-1],:) =  sizes;
width([2:2:2*m],:) =  sizes;
x = zeros(2*m,1);
x([1:2:2*m-1],:) = [0:m-1]';
x([2:2:2*m],:) = [1:m]';

% first solution
subplot(411)
hold off
plot([x;flipud(x);0], [0.5*(width(1,1)-width(:,1)); ...
    flipud(0.5*(width(1,1)+width(:,1))); 0]);
axis([-0.2, m+0.2, -0.1, 0.1+width(1,1)]);
Pos=get(gca,'Position');
pos4 = Pos(4);  % height of figure 1
pos2 = Pos(2);  % distance from bottom of subfigure to bottom of figure
set(gca,'XTickLabels',[]);
hold on;
fill([x;flipud(x);0]', [0.5*(width(1,1)-width(:,1)); ...
  flipud(0.5*(width(1,1)+width(:,1))); 0]', 0.9*ones(size([x;x;0]'))); 
caxis([-1,1]);
plot([x;flipud(x);0], [0.5*(width(1,1)-width(:,1)); ...
    flipud(0.5*(width(1,1)+width(:,1))); 0]);


% second solution
subplot(412)
hold off
plot([x;flipud(x);0], [0.5*(width(1,2)-width(:,2)); ...
    flipud(0.5*(width(1,2)+width(:,2))); 0]);
set(gca,'XTickLabels',[]);
axis([-0.2, m+0.2, -0.1, 0.1+width(1,2)]);
Pos=get(gca,'Position');
% dist is distance from top of figure to bottom of first figure
dist = pos2-Pos(2)-Pos(4);   
% pos2 is distance from bottom of subfigure to bottom of figure
pos2 = Pos(2);               
% use same scale as 1st figure 
Pos(4) = pos4*(0.2+width(1,2))/(0.2+width(1,1));  
set(gca,'Position',Pos);
hold on;
fill([x;flipud(x);0]', [0.5*(width(1,2)-width(:,2)); ...
   flipud(0.5*(width(1,2)+width(:,2))); 0]', 0.9*ones(size([x;x;0]'))); 
caxis([-1,1]);
plot([x;flipud(x);0], [0.5*(width(1,2)-width(:,2)); ...
    flipud(0.5*(width(1,2)+width(:,2))); 0]);

% third solution
subplot(413)
hold off
plot([x;flipud(x);0], [0.5*(width(1,3)-width(:,3)); ...
    flipud(0.5*(width(1,3)+width(:,3))); 0]);
set(gca,'XTickLabels',[]);
axis([-0.2, m+0.2, -0.1, 0.1+width(1,3)]);
Pos=get(gca,'Position');
% use same scale as 1st figure
Pos(4) = pos4*(0.2+width(1,3))/(0.2+width(1,1));
% shift upward s.t. distance to second plot is equal to dist
Pos(2) = pos2 - dist - Pos(4);
% pos2 is distance from bottom of plot to bottom of figure
pos2=Pos(2);
set(gca,'Position',Pos);
hold on;
fill([x;flipud(x);0]', [0.5*(width(1,3)-width(:,3)); ...
  flipud(0.5*(width(1,3)+width(:,3))); 0]', 0.9*ones(size([x;x;0]'))); 
caxis([-1,1]);
plot([x;flipud(x);0], [0.5*(width(1,3)-width(:,3)); ...
    flipud(0.5*(width(1,3)+width(:,3))); 0]);

% fourth solution
subplot(414)
hold off
plot([x;flipud(x);0], [0.5*(width(1,4)-width(:,4)); ...
    flipud(0.5*(width(1,4)+width(:,4))); 0]);
axis([-0.2, m+0.2, -0.1, 0.1+width(1,4)]);
Pos=get(gca,'Position');
Pos(4) = pos4*(0.2+width(1,4))/(0.2+width(1,1));
Pos(2) = pos2 - dist - Pos(4);
set(gca,'Position',Pos);
hold on;
fill([x;flipud(x);0]', [0.5*(width(1,4)-width(:,4)); ...
   flipud(0.5*(width(1,4)+width(:,4))); 0]', 0.9*ones(size([x;x;0]'))); 
caxis([-1,1]);
plot([x;flipud(x);0], [0.5*(width(1,4)-width(:,4)); ...
    flipud(0.5*(width(1,4)+width(:,4))); 0]);



%%%%%%%%%%  plot step responses for first solution  %%%%%%%%%%

figure(3)   

x = sizes(:,1); 
C = reshape(CC*[1;x],n,n);  
G = reshape(GG*[1;x],n,n);

% linear system:  C*vdot = -G*v + G*e*u;  y = CCC*v + D*u
A = -inv(C)*G;   
B = inv(C)*G*ones(n,1);
CCC = eye(n);
D = zeros(n,1);

% compute and plot step response
T = linspace(0,2000,1000);
[Y,X] = step(A,B,CCC,D,1,T);
hold off; plot(T,Y,'-');  hold on;
xlabel('time');
ylabel('v');

% compute threshold delay, elmore delay, dominant time constant
ind = find(Y(:,n)>0.5);
tthres=T(ind(1));
tdom=max(eig(inv(G)*C));
telm=max(sum((inv(G)*C)')); 
plot(tdom*[1;1], [0;1], '--', telm*[1;1], [0;1],'--', ...
     tthres*[1;1], [0;1], '--'); 
text(tdom,0,'d');
text(telm,0,'e');
text(tthres,0,'t');



%%%%%%%%%%  plot step responses for second solution  %%%%%%%%%%

figure(4)   

x = sizes(:,4);
C = reshape(CC*[1;x],n,n);  
G = reshape(GG*[1;x],n,n);

% linear system:  C*vdot = -G*v + G*e*u;  y = CCC*v + D*u
A = -inv(C)*G;   
B = inv(C)*G*ones(n,1);
CCC = eye(n);
D = zeros(n,1);

% compute and plot step response
T = linspace(0,2000,1000);
[Y,X] = step(A,B,CCC,D,1,T);
hold off; plot(T,Y,'-');  hold on;
xlabel('time');
ylabel('v');

% compute threshold delay, elmore delay, dominant time constant
ind = find((Y(:,n)>0.5));
tthres=T(ind(1));
tdom=max(eig(inv(G)*C));
telm=max(sum((inv(G)*C)')); 
plot(tdom*[1;1], [0;1], '--', telm*[1;1], [0;1],'--', ...
     tthres*[1;1], [0;1], '--'); 
text(tdom,0,'d');
text(telm,0,'e');
text(tthres,0,'t');

