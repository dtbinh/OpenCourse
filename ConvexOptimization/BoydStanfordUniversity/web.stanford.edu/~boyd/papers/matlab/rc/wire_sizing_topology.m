% WIRE_DRIVER_TOPOLOGY    Wire sizing and topology design
%                       (figures are generated)
% Section 5.3,  L. Vandenberghe, S. Boyd, and A. El Gamal  "Optimizing 
%                                   dominant time constant in RC circuits"
% Original by Lieven Vandenberghe 
% Adapted for CVX by Joelle Skaf - 11/25/05
%
% We size the wires for an interconnect circuit with four nodes. The 
% topology of the circuit is more complex; the wires don't even form a tree
% (refer to Figure 13 in the paper).
% The problem can be formulated with the following SDP:
%               minimize        sum(xi*li)
%                   s.t.        0 <= xi <= wmax
%                               Tmax*G(x) - C(x) >= 0
% Please refer to the paper (section 2) to find what G(x) and C(x) are. 

cvx_quiet(1);
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

%%%%%%%%%%  compute tradeoff curve  %%%%%%%%%%

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
    cvx_begin
        variable x(6)
        minimize (sum(x))
        x >= 0
        x <= wmax 
        delay*reshape(GG*[1;x],n,n) - reshape(CC*[1;x],n,n) == semidefinite(n)
    cvx_end
    areas(i) = sum(x);
end;

figure(1);
ind = finite(areas);
plot(areas(ind), delays(ind));
xlabel('area');
ylabel('Tdom');
title('Tradeoff curve');

%%%%%%%%%%  compute three solutions  %%%%%%%%%%

delays = [200 400 600];
sizes = zeros(m,3);

for i=[1,2,3]
    delay = delays(i);
    disp([' ']);
    disp(['Compute solution ', int2str(i), ...
        ' (delay = ', num2str(delay), ').']);
        cvx_begin
        variable x(6)
        minimize (sum(x))
        x >= 0
        x <= wmax 
        delay*reshape(GG*[1;x],n,n) - reshape(CC*[1;x],n,n) == semidefinite(n)
    cvx_end
    sizes(:,i) = x;
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
