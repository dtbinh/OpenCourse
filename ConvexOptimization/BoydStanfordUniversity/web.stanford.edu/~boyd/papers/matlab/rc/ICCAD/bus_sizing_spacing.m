% BUS_SIZING_SPACING        Sizing of Buses with Crosstalk 
%                           (figures are generated)
% Section 5, L. Vandenberghe, S. Boyd, and A. El Gamal  "Optimal Wire and
%                    Transistor Sizing for Circuits with Non-Tree Topology" 
% Original by Lieven Vandenberghe
% Adapted to CVX by Argyris Zymnis - 12/10/05
% 
% We consider the problem of determining the optimal line widths and
% spacings for a bus taking into account the coupling capacitances 
% between the lines. We consider an example with three wires, each
% consisting of 5 segments. The problem parameters are the segment
% widths for each wire (wij), the spacing between the three
% wires (s1, s2) and the inverse of the spacings between individual
% segments (tij). The problem can be
% formulated as the following SDP:
%               minimize        s1+s2
%                   s.t.        Tmax*G(w) - C(w,t) >= 0
%                               1/t1j <= s1 - w1j - 0.5*w2j
%                               1/t2j <= s2 - w3j - 0.5*w2j
%                               wmin <= wij <= wmax
%                               0 <= 1/tij <= smin
%
% The second and third constraints can be formulated as LMIs


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




% two points on the tradeoff curve
delays = [130 90];
sizes = [];

for i = 1:2
    fprintf(1,'Solving instance %d of 2...\n',i);
    delay = delays(i);
    cvx_begin
        variables w(3*m) t(2*m) s1(1) s2(1)
        minimize(s1+s2)
        subject to
            delay*reshape(GG*[1;w;zeros(2*m+2,1)],N,N)-reshape(CC*[1;w;t;0;0],N,N) == semidefinite(N)
            
            w >= wmin
            w <= wmax
            t >= 0
            t <= 1/smin
            
            % LMIS for inequalities
            for i = 1:m
                LMI1 = [t(i)+s1-w(i)-0.5*w(m+i) 0                       2                      ;
                        0                       t(i)+s1-w(i)-0.5*w(m+i) t(i)-s1+w(i)+0.5*w(m+i);
                        2                       t(i)-s1+w(i)+0.5*w(m+i) t(i)+s1-w(i)-0.5*w(m+i)];
                LMI1 == semidefinite(3)
                LMI2 = [t(m+i)+s2-w(2*m+i)-0.5*w(m+i) 0                             2                            ;
                        0                             t(m+i)+s2-w(2*m+i)-0.5*w(m+i) t(m+i)-s2+w(2*m+i)+0.5*w(m+i);
                        2                             t(m+i)-s2+w(2*m+i)+0.5*w(m+i) t(m+i)+s2-w(2*m+i)-0.5*w(m+i)];
                LMI2 == semidefinite(3)
            end
    cvx_end
    sizes = [sizes [w;t;s1;s2]];
end


%%%%%%%%%%  show first solution  %%%%%%%%%%

figure(1); 
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

set(gca,'XTick',[0 1 2 3 4 5 ]);
set(gca,'YTick',[0 1 2 3 4 5 6]);
title('Solution 1.');


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
figure(2); hold off;

T = linspace(0,2*delays(1),1000);
[Y1,X1] = step(A,B,CCC,D,1,T);
hold off; plot(T,Y1(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y1(1000,n),'v1');
text(T(1000),Y1(1000,2*n),'v2');
text(T(1000),Y1(1000,3*n),'v3');
axis([0 2*delays(1) -0.1 1.1]); 

% show dominant time constant
plot(delays(1)*[1;1], [-0.1;1.1], '--');

set(gca,'XTick',[0 100 200]);
set(gca,'YTick',[0 0.5 1]);
title('1st solution.  Step at 1st input.');



% response to step at 2nd input
figure(3);  hold off;

[Y2,X2] = step(A,B,CCC,D,2,T);
hold off; plot(T,Y2(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y2(1000,n),'v1');
text(T(1000),Y2(1000,2*n),'v2');
text(T(1000),Y2(1000,3*n),'v3');
axis([0 2*delays(1) -0.1 1.1]); 

plot(delays(1)*[1;1], [-0.1;1.1], '--');

set(gca,'XTick',[0 100 200]);
set(gca,'YTick',[0 0.5 1]);
title('1st solution.  Step at 2nd input.');


% response to step at 3rd input
figure(4);  hold off;

[Y3,X3] = step(A,B,CCC,D,3,T);
hold off; plot(T,Y3(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y3(1000,n),'v1');
text(T(1000),Y3(1000,2*n),'v2');
text(T(1000),Y3(1000,3*n),'v3');
axis([0 2*delays(1) -0.1 1.1]); 

plot(delays(1)*[1;1], [-0.1;1.1], '--');

set(gca,'XTick',[0 100 200]);
set(gca,'YTick',[0 0.5 1]);
title('1st solution.  Step at 3rd input.');



%%%%%%%%%%  show second solution  %%%%%%%%%%

figure(5); hold off;
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

set(gca,'XTick',[0 1 2 3 4 5 ]);
set(gca,'YTick',[0 2 4 6 8 10 12 14]);
title('Solution 2.');


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
figure(6); hold off;

T = linspace(0,2*delays(2),1000);
[Y1,X1] = step(A,B,CCC,D,1,T);
hold off; plot(T,Y1(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y1(1000,n),'v1');
text(T(1000),Y1(1000,2*n),'v2');
text(T(1000),Y1(1000,3*n),'v3');
axis([0 2*delays(2) -0.1 1.1]); 

% show dominant time constant
plot(delays(2)*[1;1], [-0.1;1.1], '--');

set(gca,'XTick',[0 50 100 150]);
set(gca,'YTick',[0 0.5 1]);
title('2nd solution.  Step at 1st input.');




% response to step at 2nd input
figure(7);  hold off;

[Y2,X2] = step(A,B,CCC,D,2,T);
hold off; plot(T,Y2(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y2(1000,n),'v1');
text(T(1000),Y2(1000,2*n),'v2');
text(T(1000),Y2(1000,3*n),'v3');
axis([0 2*delays(2) -0.1 1.1]); 

plot(delays(2)*[1;1], [-0.1;1.1], '--');
text(delays(2),-0.1,'T');
set(gca,'XTick',[0 50 100 150]);
set(gca,'YTick',[0 0.5 1]);
title('2nd solution.  Step at 2nd input.');



% response to step at 3rd input
figure(8);  hold off;

[Y3,X3] = step(A,B,CCC,D,3,T);
hold off; plot(T,Y3(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y3(1000,n),'v1');
text(T(1000),Y3(1000,2*n),'v2');
text(T(1000),Y3(1000,3*n),'v3');
axis([0 2*delays(2) -0.1 1.1]); 

plot(delays(2)*[1;1], [-0.1;1.1], '--');

set(gca,'XTick',[0 50 100 150]);
set(gca,'YTick',[0 0.5 1]);
title('2nd solution.  Step at 3rd input.');
