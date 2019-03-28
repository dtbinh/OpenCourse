% WIRE_DRIVER_SIZING    Combined sizing of drivers, repeaters, and wire
%                       (figures are generated)
% Section 5.2,  L. Vandenberghe, S. Boyd, and A. El Gamal  "Optimizing 
%                                   dominant time constant in RC circuits"
% Original by Lieven Vandenberghe 
% Adapted for CVX by Joelle Skaf - 11/25/05
%
% The first driver drives an interconnect wire, modeled as n RC Pi segments
% connected to a repeater, which drives a capacitive load through another n
% segment wires. The problem is to determine the sizes of the wire segments
% (x1, . . . , x40) and the sizes of the driver & repeater d1 and d2.
% We want to minimize area subject to bound on the combined delay Tdom1 + 
% Tdom2 of the two stages.
%               minimize        L(d1 + d2) + sum(xi*li)
%                   s.t.        0 <= xi <= wmax
%                               d1 >=0 , d2 >= 0
%                               (Tmax/2)G1(x, d1, d2) - C1(x,d2) >= 0
%                               (Tmax/2)G2(x, d1, d2) - C2(x) >= 0



cvx_quiet(1);
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


%%%%%%%%%%  compute tradeoff curve  %%%%%%%%%%

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
    cvx_begin
        variables x1(m) x2(m) d1(1) d2(1) 
        minimize ( L*(d1 + d2) + sum(x1) + sum(x2))
        x1 >= 0 
        x2 >= 0 
        d1 >= 0
        d2 >= 0
        x1 <= wmax
        x2 <= wmax
        (delay/2)*reshape(GG1*[1;x1;x2;d1;d2],n,n) - reshape(CC1*[1;x1;x2;d1;d2],n,n) == semidefinite(n)
        (delay/2)*reshape(GG2*[1;x1;x2;d1;d2],n,n) - reshape(CC2*[1;x1;x2;d1;d2],n,n) == semidefinite(n)
    cvx_end
    % optimal value
    areas(i) = cvx_optval;
end

figure(1);
ind = finite(areas);
plot(areas(ind), delays(ind));
xlabel('area');
ylabel('Tdom');
title('Area-delay tradeoff');

%%%%%%%%%%  draw one solution  %%%%%%%%%% 

delay = 189;
disp([' ']);
disp(['Compute solution for delay = ', num2str(delay), '.']);
cvx_begin
    variables x1(m) x2(m) d1(1) d2(1)
    minimize ( L*(d1 + d2) + sum(x1) + sum(x2))
    x1 >= 0
    x2 >= 0
    d1 >= 0
    d2 >= 0
    x1 <= wmax
    x2 <= wmax
    (delay/2)*reshape(GG1*[1;x1;x2;d1;d2],n,n) - reshape(CC1*[1;x1;x2;d1;d2],n,n) == semidefinite(n)
    (delay/2)*reshape(GG2*[1;x1;x2;d1;d2],n,n) - reshape(CC2*[1;x1;x2;d1;d2],n,n) == semidefinite(n)
cvx_end

figure(2)   
colormap(gray)

os = 3;  % plot drivers as block with width os and height L/(2*os);

% contour of  first wire
width1 = zeros(2*m,1);
width1([1:2:2*m-1],:) = x1;
width1([2:2:2*m],:) = x1;
x_1 = zeros(2*m,1);
x_1([1:2:2*m-1],:) = [0:m-1]';
x_1([2:2:2*m],:) = [1:m]';

% contour of second wire
width2 = zeros(2*m,1);
width2([1:2:2*m-1],:) = x2;
width2([2:2:2*m],:) = x2;
x_2 = zeros(2*m,1);
x_2([1:2:2*m-1],:) = [m:2*m-1]';
x_2([2:2:2*m],:) = [m+1:2*m]';

% drivers
xd1=[0;os;os;0;0];  
xd2=[os+m+0.5; os+m+0.5+os; os+m+0.5+os; os+m+0.5; os+m+0.5];


% draw four contours with 0.5 space between them
hold off
plot((os+0.5)+[x_1;flipud(x_1);0], ...
     [-0.5*width1; flipud(0.5*width1); -0.5*width1(1)], '-',  ...
     (os+1.5+os)+[x_2;flipud(x_2);m], ...
     [-0.5*width2; flipud(0.5*width2); -0.5*width2(1)],'-', ...
     xd1, (L/(2*os))*[-d1;-d1;d1;d1;-d1], '-', ...
     0.5+xd2, (L/(2*os))*[-d2;-d2;d2;d2;-d2], '-'); 
hold on;

% shade
fill((os+0.5)+[x_1;flipud(x_1);0], ...
     [-0.5*width1; flipud(0.5*width1); -0.5*width1(1)], ... 
     0.9*ones(size([x_1;x_1;0]')));
fill((os+1.5+os)+[x_2;flipud(x_2);m], ...
     [-0.5*width2; flipud(0.5*width2); -0.5*width2(1)], ...
     0.9*ones(size([x_2;x_2;0]')));
fill(xd1, (L/(2*os))*[-d1;-d1;d1;d1;-d1], 0.8*ones(size(xd1)));
fill(xd2+0.5,(L/(2*os))*[-d2;-d2;d2;d2;-d2], 0.8*ones(size(xd2))); 
plot((os+0.5)+[x_1;flipud(x_1);0], ...
         [-0.5*width1; flipud(0.5*width1); -0.5*width1(1)], '-',  ...
     (os+1.5+os)+[x_2;flipud(x_2);m], ...
         [-0.5*width2; flipud(0.5*width2); -0.5*width2(1)],'-', ...
     xd1, (L/(2*os))*[-d1;-d1;d1;d1;-d1], '-', ...
     xd2+0.5, (L/(2*os))*[-d2;-d2;d2;d2;-d2], '-'); 
caxis([-1,1]);
title('Solution for the point on the tradeoff curve');

% reset axes
axis([-0.1 os+os+2*m+1.6 -5 5]);
set(gca,'XTick',[os+0.5 os+ceil(m/2)+.5 os+m+0.5 ...
    2*os+m+1.5 2*os+m+ceil(m/2)+1.5 2*os+2*m+1.5]);
set(gca,'XTickLabel',[' 0'; sprintf('%2d',ceil(m/2)); ...
    sprintf('%2d',m); ' 0'; sprintf('%2d',ceil(m/2)); ...
    sprintf('%2d',m)]);

%%%%%%%%%%  step responses  %%%%%%%%%%

figure(3); 
hold off;

% first circuit
C1 = reshape(CC1*[1;x1;x2;d1;d2],n,n);
G1 = reshape(GG1*[1;x1;x2;d1;d2],n,n);

% linear system:  C1*vdot = -G1*v + G1*e*u;   y = v 
A1 = -inv(C1)*G1;   
B1 = inv(C1)*G1*ones(n,1);
CCC1 = eye(n);
D1 = zeros(n,1);

% dominant time constant and Elmore delay first wire
tdom1=max(eig(inv(G1)*C1));
telm1=max(sum((inv(G1)*C1)'));

% second circuit
C2 = reshape(CC2*[1;x1;x2;d1;d2],n,n);
G2 = reshape(GG2*[1;x1;x2;d1;d2],n,n);

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
title('Step responses for the solution');

