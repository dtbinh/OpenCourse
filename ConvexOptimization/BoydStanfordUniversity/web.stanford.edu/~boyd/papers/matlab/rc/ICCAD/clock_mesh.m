% CLOCK_MESH        Sizing of Clock Meshes 
%                   (figures are generated)
% Section 4, L. Vandenberghe, S. Boyd, and A. El Gamal  "Optimal Wire and
%                    Transistor Sizing for Circuits with Non-Tree Topology" 
% Original by Lieven Vanderberghe
% Adapted to CVX by Argyris Zymnis - 12/04/05
% 
% We consider the problem of sizing a clock mesh, so as to minimize the
% total dissipated power under a constraint on the dominant time constant.
% The numbers of nodes in the mesh is N per row or column (thus n=(N+1)^2
% in total). We divide the wire into m segments of width xi, i = 1,...,m 
% which is constrained as 0 <= xi <= Wmax. We use a pi-model of each wire 
% segment, with capacitance beta_i*xi and conductance alpha_i*xi.
% Defining C(x) = C0+x1*C1+x2*C2+...+xm*Cm we have that the dissipated 
% power is equal to ones(1,n)*C(x)*ones(n,1). Thus to minimize the
% dissipated power subject to a constraint in the widths and a constraint
% in the dominant time constant, we solve the SDP
%               minimize        ones(1,m)*C(x)*ones(m,1)
%                   s.t.        Tmax*G(x) - C(x) >= 0
%                               0 <= xi <= Wmax

clear
cvx_quiet(true);

% circuit parameters
dim=4;           % grid is dimxdim (assume dim is even)
n=(dim+1)^2;     % number of nodes
m=2*dim*(dim+1); % number of wires 
                 % 1...dim(dim+1) are horizontal segments 
                 % (numbered rowwise);
                 % dim(dim+1)+1 ... 2*dim(dim+1) are vertical
                 % (numbered columnwise)
beta = 0.5;      % capacitance per segment is twice beta times xi
alpha = 1;       % conductance per segment is alpha times xi
G0 = 1;          % source conductance
C0 = reshape([    10     2     7     5     3;
                   8     3     9     5     5;
                   1     8     4     9     3;
                   7     3     6     8     2;
                   5     2     1     9    10 ], n, 1);
wmax = 1;       % upper bound on x


% capacitance matrix 
CC = zeros(n*n,m+1);

% constant term
CC(:,1) = reshape(diag(C0),n*n,1);


% capacitances from horizontal segments 
for i=1:dim+1    % loop over rows
   for j=1:dim   % loop over segments of row i
      node1=(i-1)*(dim+1)+j;    % left node
      node2=(i-1)*(dim+1)+j+1;  % right node
      CC([n*(node1-1)+node1; n*(node2-1)+node2],1+j+(i-1)*dim) ...
         = beta*[1;1];
   end;
end;


% capacitances from columns segments 
for i=1:dim+1    % loop over columns
   for j=1:dim   % loop over segments of column i
      node1=(j-1)*(dim+1)+i;    % top node   
      node2=j*(dim+1)+i;        % bottom node
      CC([n*(node1-1)+node1; n*(node2-1)+node2],1+dim*(dim+1)+(i-1)*dim+j) ...
         = beta*[1;1];
   end;
end;



% conductance matrix without source conductances
GG = zeros(n*n,m+1);

% constant term (source conductance)
for i=1:dim+1
  node = dim/2+1 + (i-1)*(dim+1);   % first driver (middle of row 1)
  GG((node-1)*n+node,1) = G0;
end;


% loop over horizontal segments
for i=1:dim+1    % loop over rows
   for j=1:dim   % loop over segments of row i
      node1=(i-1)*(dim+1)+j;    % left node
      node2=(i-1)*(dim+1)+j+1;  % right node
      GG([n*(node1-1)+node1; n*(node1-1)+node2; ...
          n*(node2-1)+node1; n*(node2-1)+node2], 1+j+(i-1)*dim) ... 
          = alpha*[1;-1;-1;1];
   end;
end;


% loop over vertical segments
for i=1:dim+1    % loop over columns
   for j=1:dim   % loop over segments of column i
      node1=(j-1)*(dim+1)+i;    % top node   
      node2=j*(dim+1)+i;        % bottom node
      GG([n*(node1-1)+node1; n*(node1-1)+node2; ...
          n*(node2-1)+node1; n*(node2-1)+node2], 1+dim*(dim+1)+(i-1)*dim+j) ...
          = alpha*[1;-1;-1;1];
   end;
end;

%%%%%%%%%%  compute tradeoff curve  %%%%%%%%%%

% points on the tradeoff curve
nopts=20; 
delays = linspace(50,150,nopts);
areas = zeros(1,nopts);


% compute tradeoff curve
for i=1:nopts

   delay = delays(i);

   disp([' ']);
   disp(['Point ', int2str(i), ...
        ' on the tradeoff curve (delay = ', num2str(delay), ').']);

   cvx_begin
        variable x(m)
        minimize(ones(1,n)*reshape(CC*[1;x],n,n)*ones(n,1))
        subject to
            delay*reshape(GG*[1;x],n,n)-reshape(CC*[1;x],n,n) == semidefinite(n)
            x>=0
            x<=wmax
   cvx_end
        

   % optimal value
   areas(i) = sum(x);


end;


figure(1);
ind = finite(areas);
plot(2*beta*areas(ind)+sum(sum(C0)), delays(ind));
xlabel('power');
ylabel('Tdom');
axis([135 150 40 160])

%%%%%%%%%%  two solutions  %%%%%%%%%%

delays_2sol = [50; 100]
sizes = zeros(m,2);

for i=[1,2]

   delay = delays_2sol(i);

   disp([' ']);
   disp(['Compute solution ', int2str(i), ...
        ' (delay = ', num2str(delay), ').']);

   cvx_begin
        variable x(m)
        minimize(ones(1,n)*reshape(CC*[1;x],n,n)*ones(n,1))
        subject to
            delay*reshape(GG*[1;x],n,n)-reshape(CC*[1;x],n,n) == semidefinite(n)
            x>=0
            x<=wmax
   cvx_end
   
   sizes(:,i) = x;

end;


disp(['Solution 1:']);
disp(['vertical segments']); %assuming clock drivers are the middle row
reshape(sizes(1:dim*(dim+1),1),dim,dim+1)
disp(['horizontal segments']);
reshape(sizes(dim*(dim+1)+[1:dim*(dim+1)],1),dim,dim+1)'

disp(['Solution 2:']);
disp(['vertical segments']); %assuming clock drivers are the middle row
reshape(sizes(1:dim*(dim+1),2),dim,dim+1)
disp(['horizontal segments']);
reshape(sizes(dim*(dim+1)+[1:dim*(dim+1)],2),dim,dim+1)'


%%%%%%%%%%  step responses  %%%%%%%%%%


% step responses for x1 
x = sizes(:,1);

figure(2)  

C = reshape(CC(1:n*n,:)*[1;x],n,n);
G = reshape(GG(1:n*n,:)*[1;x],n,n); 

A = -inv(C)*G;
B = inv(C)*G*ones(n,1);
CCC = eye(n);
D = zeros(n,1);
  
T = linspace(0,500,2000);
[Y,X] = step(A,B,CCC,D,1,T);
indmax=0; indmin=Inf;
for j=1:size(Y,2);  
   inds = find(Y(:,j) >= 0.5);
   if (inds(1)>indmax)
      indmax = inds(1);
      jmax=j;
   end;
   if (inds(1)<indmin)
      indmin = inds(1);
      jmin=j;
   end;
end;
tthres = T(indmax);
tdom=max(eig(inv(G)*C));
elmore=max(sum((inv(G)*C)'));

hold off; plot(T,Y(:,jmax),'-',T,Y(:,jmin));  hold on;
plot(tdom*[1;1], [0;1],'--', elmore*[1;1], [0;1],'--', ...
   tthres*[1;1], [0;1], '--');
axis([0 500 0 1])
text(tdom,1,'d');
text(elmore,1,'e');
text(tthres,1,'t');
text(T(600),Y(600,jmax),['v', int2str(jmax)]);
text(T(600),Y(600,jmin),['v', int2str(jmin)]);
title('Solution 1.');



% step responses for x2 
x = sizes(:,2);

 
figure(3);
C = reshape(CC(1:n*n,:)*[1;x],n,n);
G = reshape(GG(1:n*n,:)*[1;x],n,n); 

A = -inv(C)*G;
B = inv(C)*G*ones(n,1);
CCC = eye(n);
D = zeros(n,1);
  
T = linspace(0,500,2000);
[Y,X] = step(A,B,CCC,D,1,T);
indmax=0; indmin=Inf;
for j=1:size(Y,2)  
   inds = find(Y(:,j) >= 0.5);
   if (inds(1)>indmax)
      indmax = inds(1);
      jmax=j;
   end;
   if (inds(1)<indmin)
      indmin = inds(1);
      jmin=j;
   end;
end;
tthres = T(indmax);
tdom=max(eig(inv(G)*C));
elmore=max(sum((inv(G)*C)'));

hold off; plot(T,Y(:,jmax),'-',T,Y(:,jmin));  hold on;
axis([0 500 0 1])
plot(tdom*[1;1], [0;1],'--', elmore*[1;1], [0;1],'--', ...
   tthres*[1;1], [0;1], '--');
text(tdom,1,'d');
text(elmore,1,'e');
text(tthres,1,'t');
text(T(600),Y(600,jmax),['v', int2str(jmax)]);
text(T(600),Y(600,jmin),['v', int2str(jmin)]);
title('Solution 2.');

