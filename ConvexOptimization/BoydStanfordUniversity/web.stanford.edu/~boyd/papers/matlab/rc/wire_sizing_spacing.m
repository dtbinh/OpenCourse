% WIRE_SIZING_SPACING       Combined wire sizing and spacing
%                           (figures are generated)
% Section 5.5,  L. Vandenberghe, S. Boyd, and A. El Gamal  "Optimizing 
%                                   dominant time constant in RC circuits"
% Original by Lieven Vandenberghe 
% Adapted for CVX by Joelle Skaf - 11/27/05
%
% The problem is to determine the optimal sizes of interconnect wires and 
% the optimal distances between them. We will consider an example with 3 
% wires, each consisting of 5 segments (see paper, fig.21). The variables 
% are the widths wij , and the distances s1 and s2 between the wires.
% The difference with the models used in other scripts is that we include a
% parasitic capacitance between the wires. 
% The objective is to minimize the total width s1+s2.
% The problem can be formulated with the following SDP:
%               mimimize    s1 + s2
%                    s.t.   Tmax*G(w11,..,w35)-C(w11,..,w35,t11,..,t23) >=0
%                           1/t1j <= s1 - w1j - 0.5*w2j     , j = 1,..,5
%                           1/t2j <= s2 - w3j - 0.5*w2j     , j = 1,..,5
%                           0 <= tij <= 1         , i = 1,2 , j = 1,..,5
%                           t1 >=0, t2 >= 0, s1 >=smin, s2>=smin 
%                           0 <= wij <= wmax
% the 2nd and 3rd constraints are nonlinear convex constraints that can be 
% cast as 3 x 3-LMIs (Please refer to the paper for more details).

cvx_quiet(1);
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


%%%%%%%%%%  compute two solutions  %%%%%%%%%%

% points on tradeoff curve
nopts = 50;
delays = linspace(85,200,nopts);
areas = zeros(1,nopts);

% compute tradeoff curve
for j=1:nopts
    delay = delays(j);
    disp([' ']);
    disp(['Point ', int2str(j), ...
        ' on the tradeoff curve (delay = ', num2str(delay), ').']);

    cvx_begin
        variables w1(m) w2(m) w3(m) t1(m) t2(m) s1(1) s2(1)
        minimize (s1 + s2)
        delay*reshape(GG*[1;w1;w2;w3;t1;t2;s1;s2],N,N) - reshape(CC*[1;w1;w2;w3;t1;t2;s1;s2],N,N) == semidefinite(N)
        w1 >= wmin
        w2 >= wmin
        w3 >= wmin
        t1 >= 0
        t2 >= 0
        s1 >= 0
        s2 >= 0
        w1 <= wmax
        w2 <= wmax
        w3 <= wmax
        t1 <= 1/smin
        t2 <= 1/smin
        s1 <= smax
        s2 <= smax
        for i=1:m
            LMI1 = [ t1(i)+s1-w1(i)-0.5*w2(i)  0                         2                      ;...
                0                         t1(i)+s1-w1(i)-0.5*w2(i) t1(i)-s1+w1(i)+0.5*w2(i);...
                2                         t1(i)-s1+w1(i)+0.5*w2(i) t1(i)+s1-w1(i)-0.5*w2(i)];
            LMI1 == semidefinite(3)
            LMI2 = [ t2(i)+s2-w3(i)-0.5*w2(i)  0                         2                      ;...
                0                         t2(i)+s2-w3(i)-0.5*w2(i) t2(i)-s2+w3(i)+0.5*w2(i);...
                2                         t2(i)-s2+w3(i)+0.5*w2(i) t2(i)+s2-w3(i)-0.5*w2(i)];
            LMI2 == semidefinite(3)
        end
    cvx_end
    % optimal value
    areas(j) = s1+s2;
end

figure(1);
ind = finite(areas);
plot(areas(ind), delays(ind));
xlabel('total width s_1 + s_2');
ylabel('dominant time constant');
title('Tradeoff curve')

%%%%%%%%%%  compute two solutions  %%%%%%%%%%

% two points on the tradeoff curve
delays = [127 87];
sizes = zeros(5*m+2,2);

for j=[1,2]
    delay = delays(j);
    disp([' ']);
    disp(['Compute solution ', int2str(j), ...
        ' (delay = ', num2str(delay), ').']);
    cvx_begin
        variables w1(m) w2(m) w3(m) t1(m) t2(m) s1(1) s2(1)
        minimize (s1 + s2)
        delay*reshape(GG*[1;w1;w2;w3;t1;t2;s1;s2],N,N) - reshape(CC*[1;w1;w2;w3;t1;t2;s1;s2],N,N) == semidefinite(N)
        w1 >= wmin
        w2 >= wmin
        w3 >= wmin
        t1 >= 0
        t2 >= 0
        s1 >= 0
        s2 >= 0
        w1 <= wmax
        w2 <= wmax
        w3 <= wmax
        t1 <= 1/smin
        t2 <= 1/smin
        s1 <= smax
        s2 <= smax
        for i=1:m
            LMI1 = [ t1(i)+s1-w1(i)-0.5*w2(i)  0                         2                      ;...
                0                         t1(i)+s1-w1(i)-0.5*w2(i) t1(i)-s1+w1(i)+0.5*w2(i);...
                2                         t1(i)-s1+w1(i)+0.5*w2(i) t1(i)+s1-w1(i)-0.5*w2(i)];
            LMI1 == semidefinite(3)
            LMI2 = [ t2(i)+s2-w3(i)-0.5*w2(i)  0                         2                      ;...
                0                         t2(i)+s2-w3(i)-0.5*w2(i) t2(i)-s2+w3(i)+0.5*w2(i);...
                2                         t2(i)-s2+w3(i)+0.5*w2(i) t2(i)+s2-w3(i)-0.5*w2(i)];
            LMI2 == semidefinite(3)
        end
    cvx_end
    sizes(:,j) = [w1;w2;w3;t1;t2;s1;s2];
end

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
title('Solution marked (a) on the tradeoff curve');

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
title('Voltage at the output nodes due to a step at 1^{st} wire');

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
title('Voltage at the output nodes due to a step at 2^{nd} wire');

% response to step at 3rd input
figure(5);  hold off;

[Y3,X3] = step(A,B,CCC,D,3,T);
hold off; plot(T,Y3(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y3(1000,n),'v1');
text(T(1000),Y3(1000,2*n),'v2');
text(T(1000),Y3(1000,3*n),'v3');
axis([0 2*delays(1) -0.1 1.1]); 
title('Voltage at the output nodes due to a step at 3^{rd} wire');
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
title('Voltage at the output nodes due to a step at 1^{st} & 2^{nd} wires');

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
title('Voltage at the output nodes due to a step at 1^{st} & 3^{rd} wires');

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
title('Voltage at the output nodes due to a step at 2^{nd} & 3^{rd} wires');

%%%%%%%%%%  show second solution  %%%%%%%%%%

figure(9); hold off;
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
title('Solution marked (b) on the tradeoff curve');



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
figure(10); hold off;

T = linspace(0,2*delays(2),1000);
[Y1,X1] = step(A,B,CCC,D,1,T);
hold off; plot(T,Y1(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y1(1000,n),'v1');
text(T(1000),Y1(1000,2*n),'v2');
text(T(1000),Y1(1000,3*n),'v3');
axis([0 2*delays(2) -0.1 1.1]); 
title('Voltage at the output nodes due to a step at 1^{st} wire');
% show dominant time constant
plot(delays(2)*[1;1], [-0.1;1.1], '--');


% response to step at 2nd input
figure(11);  hold off;

[Y2,X2] = step(A,B,CCC,D,2,T);
hold off; plot(T,Y2(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y2(1000,n),'v1');
text(T(1000),Y2(1000,2*n),'v2');
text(T(1000),Y2(1000,3*n),'v3');
axis([0 2*delays(2) -0.1 1.1]); 
title('Voltage at the output nodes due to a step at 2^{nd} wire');
plot(delays(2)*[1;1], [-0.1;1.1], '--');
text(delays(2),-0.1,'T');


% response to step at 3rd input
figure(12);  hold off;

[Y3,X3] = step(A,B,CCC,D,3,T);
hold off; plot(T,Y3(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y3(1000,n),'v1');
text(T(1000),Y3(1000,2*n),'v2');
text(T(1000),Y3(1000,3*n),'v3');
axis([0 2*delays(2) -0.1 1.1]); 
title('Voltage at the output nodes due to a step at 3^{rd} wire');
plot(delays(2)*[1;1], [-0.1;1.1], '--');


% response to step at 1st and 2nd input
figure(13); hold off;

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
title('Voltage at the output nodes due to a step at 1^{st} & 2^{nd} wires');

% response to step at 1st and 3rd input
figure(14);  hold off;

B = inv(C)*G*[ones(n,1); zeros(n,1); ones(n,1)];

[Y5,X5] = step(A,B,CCC,D,1,T);
hold off; plot(T,Y5(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y5(1000,n),'v1');
text(T(1000),Y5(1000,2*n),'v2');
text(T(1000),Y5(1000,3*n),'v3');
axis([0 2*delays(2) -0.1 1.1]); 
plot(delays(2)*[1;1], [-0.1;1.1], '--');
title('Voltage at the output nodes due to a step at 1^{st} & 3^{rd} wires');

% response to step at 2nd and 3rd wire 
figure(15);  hold off;
B = inv(C)*G*[zeros(n,1); ones(2*n,1)];

[Y6,X6] = step(A,B,CCC,D,1,T);
hold off; plot(T,Y6(:,[n,2*n,3*n]),'-');  hold on;
text(T(1000),Y6(1000,n),'v1');
text(T(1000),Y6(1000,2*n),'v2');
text(T(1000),Y6(1000,3*n),'v3');
axis([0 2*delays(2) -0.1 1.1]); 
plot(delays(2)*[1;1], [-0.1;1.1], '--');
title('Voltage at the output nodes due to a step at 2^{nd} & 3^{rd} wires');
