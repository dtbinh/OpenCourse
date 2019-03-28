% 
% FIRST EXAMPLE in ICCAD paper:  Clock mesh  with nonsymmetric loads
%


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



% LMI constraints
%  (-1/delay)*C + G (+tI) >= 0, i=1,...,m   (nxn)
%  x >= 0
%  x <= wmax

FF = zeros(n*n+2*m, m+2);
blck_szs = [n; ones(2*m,1)];
FF(1:n*n,m+2) = reshape(eye(n),n*n,1);
FF(n*n+[1:m],1+[1:m]) = eye(m);
FF(n*n+m+[1:m],1+[1:m]) = -eye(m);


% objective: minimize power
c=ones(m,1); 




%%%%%%%%%%  compute tradeoff curve  %%%%%%%%%%

% tolerances for SP
abstol = 1e-5;
reltol = 1e-6;

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

   % LMI
   FF(1:n*n,1:m+1) = GG - (1/delay)*CC;
   
   % try very thick wires
   x0 = 0.99*wmax*ones(m,1); 
   t = min(eig(reshape(FF(1:n*n,1:m+1)*[1;x0],n,n)));
   feasible = 1;

   if (t <= 0)   % x0 is infeasible, run phase 1

      disp([' ']);  disp(['Phase 1.']);

      % use slightly tighter bounds on x to have strict feasibility
      FF(n*n+[1:m],1) = -1e-5*ones(m,1);
      FF(n*n+m+[1:m],1) = (1-1e-5)*wmax*ones(m,1); 

      % dual feasible solution
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

      % dual feasible point
      z = FF(1:n*n,2:m+1)'*reshape(eye(n)/n,n*n,1) - 1;
      Z0 = [ reshape(eye(n)/n,n*n,1); -min(z,0)+1e-3;  max(z,0)+1e-3];

      % call SP
      [x,Z,ul,info,time] = sp(FF(:,1:m+1), blck_szs, c, ...
          x0, Z0, 10.0, abstol, reltol, 0.0, 100);

      % optimal value
      areas(i) = c'*x;

   end;

end;


figure(1);
ind = finite(areas);
plot(2*beta*areas(ind)+sum(sum(C0)), delays(ind));
xlabel('power');
ylabel('Tdom');
axis([135 150 40 160])
set(gca,'XTick',[135 140 145 150]);
set(gca,'YTick',[50 100 150]);




%%%%%%%%%%  two solutions  %%%%%%%%%%

% tolerances for SP
abstol = 1e-5;
reltol = 1e-6;

delays_2sol = [50; 100]
sizes = zeros(m,2);

for i=[1,2]

   delay = delays_2sol(i);

   disp([' ']);
   disp(['Compute solution ', int2str(i), ...
        ' (delay = ', num2str(delay), ').']);

   FF(1:n*n,1:m+1) = GG - (1/delay)*CC;
   
   % try very thick wires
   x0 = 0.99*wmax*ones(m,1); 
   t = min(eig(reshape(FF(1:n*n,1:m+1)*[1;x0],n,n)));
   feasible = 1;

   if (t <= 0)   % x0 is infeasible, run phase 1

      disp([' ']);  disp(['Phase 1.']);

      % use slightly tighter bounds on x to have strict feasibility
      FF(n*n+[1:m],1) = -1e-5*ones(m,1);
      FF(n*n+m+[1:m],1) = (1-1e-5)*wmax*ones(m,1); 

      % dual feasible solution
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

      % dual feasible point
      z = FF(1:n*n,2:m+1)'*reshape(eye(n)/n,n*n,1) - 1;
      Z0 = [ reshape(eye(n)/n,n*n,1); -min(z,0)+1e-3;  max(z,0)+1e-3];

      % call SP
      [sizes(:,i),Z,ul,info,time] = sp(FF(:,1:m+1), blck_szs, c, ...
          x0, Z0, 10.0, abstol, reltol, 0.0, 100);

   end;

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
set(gca,'XTick',[0 100 200 300 400 500]);
set(gca,'YTick',[0 0.5 1.0]);
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
set(gca,'XTick',[0 100 200 300 400 500]);
set(gca,'YTick',[0 0.5 1.0]);
text(tdom,1,'d');
text(elmore,1,'e');
text(tthres,1,'t');
text(T(600),Y(600,jmax),['v', int2str(jmax)]);
text(T(600),Y(600,jmin),['v', int2str(jmin)]);
title('Solution 2.');

