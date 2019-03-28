% WIRE_SIZING       problem of sizing an interconnect wire 
%                   (figures are generated)
% Section 5.1, L. Vandenberghe, S. Boyd, and A. El Gamal  "Optimizing 
%                                   dominant time constant in RC circuits" 
% Original by Lieven Vandenberghe 
% Adapted for CVX by Joelle Skaf - 11/25/05
% 
% we consider the problem of sizing an interconnect wire that connects
% a voltage source and conductance G to a capacitive load C. We divide the
% wire into n segments of length li, and width xi, i = 1,...,n, which is
% constrained as 0 <= xi <= Wmax. The total area of the interconnect wire 
% is therefore sum(li*xi). We use a pi-model of each wire segment, with 
% capacitors beta_i*xi and conductance alpha_i*xi. 
% To minimize the total area subject to the width bound and a bound Tmax on 
% dominant time constant, we solve the SDP
%               minimize        sum_{i=1}^20 xi*li
%                   s.t.        Tmax*G(x) - C(x) >= 0
%                               0 <= xi <= Wmax

cvx_quiet(1);
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
    cvx_begin
        variable x(m) 
        minimize (sum(x))
        delay*reshape(GG*[1;x],n,n) - reshape(CC*[1;x],n,n) == semidefinite(n)
        x >= 0 
        x <= wmax
    cvx_end
    % optimal value
    areas(i) = sum(x);
end


figure(1)
ind = finite(areas);   % skip infeasible values
plot(areas(ind), delays(ind));
xlabel('area');
ylabel('Tdom');
title('Area-delay tradeoff curve');

%%%%%%%%%%  compute four solutions  %%%%%%%%%%

% compute solution for these four values of delay
delays = [370, 400, 600, 1800];
sizes = zeros(m,4);

for i=1:4
    delay = delays(i);
    disp([' ']);
    disp(['Compute solution ', int2str(i), ...
        ' (delay = ', num2str(delay), ').']);
    cvx_begin
        variable x(m)
        minimize (sum(x))
        delay*reshape(GG*[1;x],n,n) - reshape(CC*[1;x],n,n) == semidefinite(n)
        x >= 0
        x <= wmax
    cvx_end
    % optimal value
    sizes(:,i) = x;
end
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
set(gca,'XTickLabel',[]);
hold on;
fill([x;flipud(x);0]', [0.5*(width(1,1)-width(:,1)); ...
  flipud(0.5*(width(1,1)+width(:,1))); 0]', 0.9*ones(size([x;x;0]'))); 
caxis([-1,1]);
plot([x;flipud(x);0], [0.5*(width(1,1)-width(:,1)); ...
    flipud(0.5*(width(1,1)+width(:,1))); 0]);
title('Solution at four points on the tradeoff curve');

% second solution
subplot(412)
hold off
plot([x;flipud(x);0], [0.5*(width(1,2)-width(:,2)); ...
    flipud(0.5*(width(1,2)+width(:,2))); 0]);
set(gca,'XTickLabel',[]);
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
set(gca,'XTickLabel',[]);
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
title('Step responses at the 21 nodes for the 1^{st} solution');


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
title('Step responses at the 21 nodes for the 4^{th} solution');


