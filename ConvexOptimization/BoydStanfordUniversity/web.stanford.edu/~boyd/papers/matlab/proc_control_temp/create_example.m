% CREATE_EXAMPLE   creates the problem data for numerical example in
%                  "Processor speed control with thermal constraints"
%                  by Mutapcic, Boyd, Murali, Atienza, De Micheli, Gupta
%
% Variables:
% n      - number of processors
% m      - number of nodes in the model (places with temp. observation)
% l      - number of edges between the nodes
%
% k      - thermal conductivity vector (R^l vector)
% A      - reduced node-incidence matrix for the thermal model (m x l)
% B      - processor power distribution matrix into nodes (m x n)
% phi    - structure defining speed-power model function phi (its value, grad, Hessian)
% Tamb   - ambient temperature value
% Tother - temperature vector due to uncontrollable sources (R^m vector)
% Tmax   - maximum allowable temperature
% smin   - minimum processor speed
% smax   - maximum processor speed
%
% Created variables:
%
% G      - thermal matrix (m x n)
% Ghat   - truncated (sparse) thermal matrix (m x n)

% almirm, 9/07

% define optimization problem constants
Tamb = 25;
Tmax = 75;
smin = 1;
smax = 3;

%
% define speed-power model in terms of phi function value, gradient, and Hessian
% (we use the well-known cube speed-power model and also provide square model)
%
phi_model = 'cube';

switch (phi_model)
  case 'square'
    phi.fval = @(s) s.^2;
    phi.grad = @(s) 2*s;
    phi.hess = @(s) spdiags(2*ones(length(s),1),0,length(s),length(s));
  case 'cube'
    phi.fval = @(s) s.^3;
    phi.grad = @(s) 3*s.^2;
    phi.hess = @(s) spdiags(6*s,0,length(s),length(s));
end

%
% generate thermal data for an 10x10 square grid of processors, with
% rectangular grid of 55x75 thermal sensors (nodes in the thermal model)
%
grid_n = 10;
node_n = 75; node_m = 55;

n = grid_n*grid_n;
m = node_m*node_n;
l = 2*m + node_m + node_n - 4;

% construct sensor and processor mesh grids
skip = 5;
[X,Y] = meshgrid(1:node_n,1:node_m);

% plot the grid with sensors and processors
figure(1), clf
dw = 0.25;
plot(X(:),Y(:),'.','Color',[0 .75 0],'LineWidth',1); hold on
for pX = [15:skip:node_n-15]
  for pY = [5:skip:node_m-5]
    plot([pX-dw pX+1+dw pX+1+dw pX-dw pX-dw],[pY-dw pY-dw pY+1+dw pY+1+dw pY-dw],'r-','Linewidth',1);
  end
end
hold off
xlabel('x'), ylabel('y')
set(gca,'FontSize',13);
axis([1 node_n 1 node_m]),
set(gca,'XTick',[5.5:skip:node_n],'YTick',[5.5:skip:node_m])
set(gca,'XTickLabel',[],'YTickLabel',[])

% create the reduced node-incidence matrix A
fprintf(1,'* creating reduced node-incidence matrix A\n');
A = sparse(m,l);
k = zeros(l,1);

% connect down the grid starting from the origin node (top left corner)
edge = 1;
for j = [1:node_n]
  for i = [(j-1)*node_m+1:j*node_m-1]
    A( i,   edge ) = -1;
    A( i+1, edge ) = +1;
    k(edge) = 1.0;
    edge = edge+1;
  end   
end

% connect across the grid starting from the origin node (top left corner)
for j = [1:node_n-1]
  for i = [(j-1)*node_m+1:j*node_m]
    A( i,        edge ) = -1;
    A( i+node_m, edge ) = +1;
    k(edge) = 1.0;
    edge = edge+1;
  end
end

% connect ground node to the outside nodes
grd_nodes = [1:node_m,(node_n-1)*node_m+1:m, ...
             node_m+1:node_m:(node_n-1)*node_m, ...
             2*node_m:node_m:m-1];
for node = grd_nodes
  A(node,edge) = -1;
  k(edge) = 1.5;
  edge = edge+1;
end


% construct processor power density matrix B
fprintf(1,'* creating processor power matrix B\n');
B = sparse(m,n);

% find top left corner of the processor (each one spans adjacent 4 sensors)
% here we shift everything by one to center processors
ind = [];
for pX = [14:skip:node_n-15]
  for pY = [5:skip:node_m-5]
     ind = [ind pX*node_m+pY];
  end
end

for proc = 1:n
    %  node no.           proc no.  amt.
    B( ind(proc),          proc ) = 0.25;
    B( ind(proc)+1,        proc ) = 0.25;
    B( ind(proc)+node_m,   proc ) = 0.25;
    B( ind(proc)+node_m+1, proc ) = 0.25;
end

%
% compute the thermal matrix G (which is a fully dense matrix)
%
fprintf(1,'* computing G\n');
G = (A*diag(k)*A')\B;

% obtain Ghat (truncated G with about 7% sparsity pattern)
fprintf(1,'* computing Ghat where threshold level is 0.15\n');
Ghat = G.*(G > 0.15);
Ghat = sparse(Ghat*diag(sum(G)./sum(Ghat)));

fprintf(1,'    - nnz(G)         = %d\n',nnz(G));
fprintf(1,'    - nnz(Ghat)      = %d\n',nnz(Ghat));
fprintf(1,'    - sparsity(Ghat) = %3.2g%%\n',nnz(Ghat)/(m*n)*100);

%
% construct Tother (temp. rise due to uncontrollable sources)
%
fprintf(1,'* creating Tother\n\n');
therm_const = 1;
therm_std   = 35;

% unknown powers beyond our control
p_other = [12 14 14 14 12 12 12 8 14 10 22 23 14 14 12 14 12 12 10 10 10]';
n_other = length(p_other);
G_other = zeros(m,n_other);
pX = [2  2  2  2 10 10  2  2  2  2  36 36 74 73 72 73 74 73 73 73 73];
pY = [2 10 20 30 30 34 34 38 45 52  10  2  4 15 25 30 34 40 44 50 53];
for k = 1:n_other
  g = therm_const*exp(-((X-pX(k)).^2 + (Y-pY(k)).^2)/therm_std);
  G_other(:,k) = g(:);
end

% Tother due to unknown sources and measurement noise
figure(2), clf
Tother = G_other*p_other;
Tother_image = reshape(Tother,node_m,node_n);
imagesc(Tother_image,[0 25])
colorbar
axis([1 node_n 1 node_m]),
set(gca,'XTick',[5.5:skip:node_n],'YTick',[5.5:skip:node_m])
set(gca,'XTickLabel',[],'YTickLabel',[])

fprintf(1,'Created problem data for the example.\n');
