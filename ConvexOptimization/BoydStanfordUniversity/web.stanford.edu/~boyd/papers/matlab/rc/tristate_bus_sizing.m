% TRISTATE_BUS_SIZING       Tri-state bus sizing and topology design
%                           (figures are generated)
% Section 5.4,  L. Vandenberghe, S. Boyd, and A. El Gamal  "Optimizing 
%                                   dominant time constant in RC circuits" 
% Original by Lieven Vandenberghe 
% Adapted for CVX by Joelle Skaf - 11/27/05
%
% We optimize a tri-state bus connecting six nodes ( The model for the bus 
% is shown in the paper, fig.16). The total wire area is sum_{i>j} lij*xij 
% The bus can be driven from any node. When node i drives the bus, the ith 
% switch is closed and the others are all open. Thus we really have six 
% different circuits, each corresponding to a given node driving the bus.
% we require that the dominant time constant of each of the six drive 
% configuration circuits has dominant time constant less than Tmax.
% The problem can be formulated with the following SDP:
%               minimize        sum_{i>j}(x_ij*l_ij)
%                   s.t.        0 <= xij <= wmax
%                               Tmax*(G(x) + GE_kk) - C(x) >= 0 , 1 <=k<= 6
% The matrix E_kk is zero except for the kth diagonal element, which is 1.

cvx_quiet(1);
% circuit parameters
n=6;         % number of nodes
m=15;        % number of wires
beta = 0.5;  % capacitance per segment is twice beta times xi*li
alpha = 1;   % conductance per segment is alpha times xi/li
G0 = 1;      % source conductance
C0 = 10;     % load capacitor
wmax = 1;    % upper bound on x


% positions of six nodes 
xpos = [ 0   1   6   8  -4  -1 ; 
      0  -1   4  -2   1   4 ] ;
X11 = repmat(xpos(1,:),n,1);
X12 = repmat(xpos(1,:)',1,n);
X21 = repmat(xpos(2,:),n,1);
X22 = repmat(xpos(2,:)',1,n);
LL = abs(X11-X12) + abs(X21-X22);
L = tril(LL);
L = L(find(L>0));

% capacitance matrix 
CC = zeros(n*n,m+1);

% constant term
CC(:,1) = C0*reshape(eye(n),n*n,1);

% capacitances from segments
CC([n+2,1],2) = beta*L(1)*[1;1];
CC([2*n+3,1],3) = beta*L(2)*[1;1];   
CC([3*n+4,1],4) = beta*L(3)*[1;1];
CC([4*n+5,1],5) = beta*L(4)*[1;1];   
CC([5*n+6,1],6) = beta*L(5)*[1;1];   
CC([2*n+3,n+2],7) = beta*L(6)*[1;1];   
CC([3*n+4,n+2],8) = beta*L(7)*[1;1];
CC([4*n+5,n+2],9) = beta*L(8)*[1;1];   
CC([5*n+6,n+2],10) = beta*L(9)*[1;1];   
CC([3*n+4,2*n+3],11) = beta*L(10)*[1;1];  
CC([4*n+5,2*n+3],12) = beta*L(11)*[1;1];
CC([5*n+6,2*n+3],13) = beta*L(12)*[1;1];
CC([4*n+5,3*n+4],14) = beta*L(13)*[1;1];
CC([5*n+6,3*n+4],15) = beta*L(14)*[1;1];  
CC([5*n+6,4*n+5],16) = beta*L(15)*[1;1];  


% conductance matrix without source conductances
GG = zeros(n*n,m+1);
GG([n+2,n+1,2,1],2) = alpha*[1;-1;-1;1]/L(1);
GG([2*n+3,2*n+1,3,1],3) = alpha*[1;-1;-1;1]/L(2);
GG([3*n+4,3*n+1,4,1],4) = alpha*[1;-1;-1;1]/L(3);
GG([4*n+5,4*n+1,5,1],5) = alpha*[1;-1;-1;1]/L(4);   
GG([5*n+6,5*n+1,6,1],6) = alpha*[1;-1;-1;1]/L(5);   
GG([2*n+3,2*n+2,n+3,n+2],7) = alpha*[1;-1;-1;1]/L(6);   
GG([3*n+4,3*n+2,n+4,n+2],8) = alpha*[1;-1;-1;1]/L(7);   
GG([4*n+5,4*n+2,n+5,n+2],9) = alpha*[1;-1;-1;1]/L(8);   
GG([5*n+6,n+6,5*n+2,n+2],10) = alpha*[1;-1;-1;1]/L(9);   
GG([3*n+4,3*n+3,2*n+4,2*n+3],11) = alpha*[1;-1;-1;1]/L(10);  
GG([4*n+5,4*n+3,2*n+5,2*n+3],12) = alpha*[1;-1;-1;1]/L(11);  
GG([5*n+6,5*n+3,2*n+6,2*n+3],13) = alpha*[1;-1;-1;1]/L(12);  
GG([4*n+5,4*n+4,3*n+5,3*n+4],14) = alpha*[1;-1;-1;1]/L(13);  
GG([5*n+6,5*n+4,3*n+6,3*n+4],15) = alpha*[1;-1;-1;1]/L(14);  
GG([5*n+6,5*n+5,4*n+6,4*n+5],16) = alpha*[1;-1;-1;1]/L(15);   

%%%%%%%%%%  compute tradeoff curve  %%%%%%%%%%

% points on the tradeoff curve
nopts=50;
delays = linspace(410,2000,nopts);
areas = zeros(1,nopts);


% compute tradeoff curve
for i=1:nopts
    delay = delays(i);
    disp([' ']);
    disp(['Point ', int2str(i), ...
        ' on the tradeoff curve (delay = ', num2str(delay), ').']);

    cvx_begin
        variable x(m)
        minimize (L'*x)
        delay*(reshape(GG(1:n*n,:)*[1;x],n,n) + diag([1 0 0 0 0 0])) - reshape(CC(1:n*n,:)*[1;x],n,n) ==semidefinite(n)
        delay*(reshape(GG(1:n*n,:)*[1;x],n,n) + diag([0 1 0 0 0 0])) - reshape(CC(1:n*n,:)*[1;x],n,n) ==semidefinite(n)
        delay*(reshape(GG(1:n*n,:)*[1;x],n,n) + diag([0 0 1 0 0 0])) - reshape(CC(1:n*n,:)*[1;x],n,n) ==semidefinite(n)
        delay*(reshape(GG(1:n*n,:)*[1;x],n,n) + diag([0 0 0 1 0 0])) - reshape(CC(1:n*n,:)*[1;x],n,n) ==semidefinite(n)
        delay*(reshape(GG(1:n*n,:)*[1;x],n,n) + diag([0 0 0 0 1 0])) - reshape(CC(1:n*n,:)*[1;x],n,n) ==semidefinite(n)
        delay*(reshape(GG(1:n*n,:)*[1;x],n,n) + diag([0 0 0 0 0 1])) - reshape(CC(1:n*n,:)*[1;x],n,n) ==semidefinite(n)
        x >= 0
        x <= wmax
    cvx_end
    % optimal value
    areas(i) = cvx_optval;
end;

figure(1);
ind = finite(areas);
plot(areas(ind), delays(ind));
xlabel('area');
ylabel('Tdom');

%%%%%%%%%%  two solutions  %%%%%%%%%%

% points on the tradeoff curve
delays = [410, 2000];
sizes = zeros(m,2);

% compute tradeoff curve
for i=[1,2]
    delay = delays(i);
    disp([' ']);
    disp(['Compute solution ', int2str(i), ...
        ' (delay = ', num2str(delay), ').']);
    cvx_begin
        variable x(m)
        minimize (L'*x)
        delay*(reshape(GG(1:n*n,:)*[1;x],n,n) + diag([1 0 0 0 0 0])) - reshape(CC(1:n*n,:)*[1;x],n,n) ==semidefinite(n)
        delay*(reshape(GG(1:n*n,:)*[1;x],n,n) + diag([0 1 0 0 0 0])) - reshape(CC(1:n*n,:)*[1;x],n,n) ==semidefinite(n)
        delay*(reshape(GG(1:n*n,:)*[1;x],n,n) + diag([0 0 1 0 0 0])) - reshape(CC(1:n*n,:)*[1;x],n,n) ==semidefinite(n)
        delay*(reshape(GG(1:n*n,:)*[1;x],n,n) + diag([0 0 0 1 0 0])) - reshape(CC(1:n*n,:)*[1;x],n,n) ==semidefinite(n)
        delay*(reshape(GG(1:n*n,:)*[1;x],n,n) + diag([0 0 0 0 1 0])) - reshape(CC(1:n*n,:)*[1;x],n,n) ==semidefinite(n)
        delay*(reshape(GG(1:n*n,:)*[1;x],n,n) + diag([0 0 0 0 0 1])) - reshape(CC(1:n*n,:)*[1;x],n,n) ==semidefinite(n)
        x >= 0
        x <= wmax
    cvx_end
    % optimal value
    sizes(:,i) = x;
end;

disp(['Two solutions:']);
sizes

%%%%%%%%%%  step responses  %%%%%%%%%%


% step responses for x1 
x = sizes(:,1);

for i=1:6;
    figure(i+1)

    C = reshape(CC(1:n*n,:)*[1;x],n,n);
    G = reshape(GG(1:n*n,:)*[1;x],n,n);
    G(i,i) = G(i,i) + G0;

    A = -inv(C)*G;
    B = inv(C)*G*ones(n,1);
    CCC = eye(n);
    D = zeros(n,1);

    T = linspace(0,1000,1000);
    [Y,X] = step(A,B,CCC,D,1,T);
    hold off; plot(T,Y,'-');  hold on;

    ind=0; for j=1:size(Y,2);
        inds = find(Y(:,j) >= 0.5);
        ind = max(inds(1),ind);
    end;
    tthres = T(ind);
    tdom=max(eig(inv(G)*C));
    elmore=max(sum((inv(G)*C)'));
    plot(tdom*[1;1], [0;1],'--', elmore*[1;1], [0;1],'--', ...
        tthres*[1;1], [0;1], '--');
    text(tdom,0,'d');
    text(elmore,0,'e');
    text(tthres,0,'t');
    ylabel('Voltage');
    title(['Step response for solution 2 when switch ' num2str(i) ' is closed']);
 end;

% step responses for x2 
x = sizes(:,2);

for i=1:6;
    
    figure(i+1);

    C = reshape(CC(1:n*n,:)*[1;x],n,n);
    G = reshape(GG(1:n*n,:)*[1;x],n,n);
    G(i,i) = G(i,i) + G0;

    A = -inv(C)*G;
    B = inv(C)*G*ones(n,1);
    CCC = eye(n);
    D = zeros(n,1);

    T = linspace(0,3000,1000);
    [Y,X] = step(A,B,CCC,D,1,T);
    hold off; plot(T,Y,'-');  hold on;

    ind=0; for j=1:size(Y,2);
        inds = find(Y(:,j) >= 0.5);
        ind = max(inds(1),ind);
    end;
    tthres = T(ind);
    tdom=max(eig(inv(G)*C));
    elmore=max(sum((inv(G)*C)'));
    plot(tdom*[1;1], [0;1],'--', elmore*[1;1], [0;1],'--', ...
        tthres*[1;1], [0;1], '--');
    text(tdom,0,'d');
    text(elmore,0,'e');
    text(tthres,0,'t');
    ylabel('Voltage');
    title(['Step response for solution 2 when switch ' num2str(i) ' is closed'])
end;
