% 
% EXAMPLE 4: Tri-state bus sizing and topology design
%


% circuit parameters
n=6;         % number of nodes
m=15;        % number of wires
beta = 0.5;  % capacitance per segment is twice beta times xi*li
alpha = 1;   % conductance per segment is alpha times xi/li
G0 = 1;      % source conductance
C0 = 10;     % load capacitor
wmax = 1;    % upper bound on x


% positions of six nodes 
X = [ 0   1   6   8  -4  -1 ; 
      0  -1   4  -2   1   4 ] ;

% lengths of the 15 interconnections
L = [ norm(X(:,2) - X(:,1), 1);     %  1
      norm(X(:,3) - X(:,1), 1);     %  2
      norm(X(:,4) - X(:,1), 1);     %  3
      norm(X(:,5) - X(:,1), 1);     %  4
      norm(X(:,6) - X(:,1), 1);     %  5
      norm(X(:,3) - X(:,2), 1);     %  6
      norm(X(:,4) - X(:,2), 1);     %  7
      norm(X(:,5) - X(:,2), 1);     %  8 
      norm(X(:,6) - X(:,2), 1);     %  9
      norm(X(:,4) - X(:,3), 1);     % 10
      norm(X(:,5) - X(:,3), 1);     % 11
      norm(X(:,6) - X(:,3), 1);     % 12
      norm(X(:,5) - X(:,4), 1);     % 13
      norm(X(:,6) - X(:,4), 1);     % 14
      norm(X(:,6) - X(:,5), 1) ];   % 15


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


% LMI constraints
%  (-1/delay)*C + G + Eii (+tI) >= 0, i=1,...,m   (nxn)
%  x >= 0
%  x >= wmax

FF = zeros(n*n*n+2*m, m+2);
blck_szs = [n*ones(n,1); ones(2*m,1)];
for j=1:n
   FF((j-1)*n*n+[1:n*n],m+2) = reshape(eye(n),n*n,1);
end;
FF(n*n*n+[1:m],1+[1:m]) = eye(m);
FF(n*n*n+m+[1:m],1+[1:m]) = -eye(m);


% objective: minimize area
c=L;


%%%%%%%%%%  compute tradeoff curve  %%%%%%%%%%

% tolerances for SP
abstol = 1e-4;
reltol = 1e-4;

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

   % n LMIs
   for j=1:n
      FF((j-1)*n*n+[1:n*n],1:m+1) = GG - (1/delay)*CC;
      FF((j-1)*n*n+(j-1)*n+j,1) = G0;
   end;
   
   % try very thick wire
   x0 = 0.99*wmax*ones(m,1); 
   t = Inf;
   for j=1:n
      t = min(t, min(eig(...
          reshape(FF((j-1)*n*n+[1:n*n],1:m+1)*[1;x0],n,n))));
   end;
   feasible = 1;

   if (t <= 0)   % x0 is infeasible, run phase 1

      disp([' ']);  disp(['Phase 1.']);

      % use slightly tighter bounds on x1, x2, d1, d2 to ensure
      % strict feasibility
      FF(n*n*n+[1:m],1) = -1e-5*ones(m,1);
      FF(n*n*n+m+[1:m],1) = (1-1e-5)*wmax*ones(m,1); 

      % dual feasible solution
      Z0 = zeros(n*n*n+2*m,1);
      for j=1:n
         Z0((j-1)*n*n+[1:n*n],1) = reshape(eye(n)/(n*n),n*n,1);
      end;
      z = FF(1:n*n*n,2:m+1)' * Z0(1:n*n*n);

      Z0(n*n*n+[1:2*m]) = [-min(z,0)+1e-3; max(z,0)+1e-3 ];

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
      FF(n*n*n+[1:m],1) = zeros(m,1);
      FF(n*n*n+m+[1:m],1) = wmax*ones(m,1); 

      % dual feasible point
      Z0 = zeros(n*n*n+2*m,1);
      for j=1:n
         Z0((j-1)*n*n+[1:n*n],1) = reshape(eye(n)/(n*n),n*n,1);
      end;
      z = FF(1:n*n*n,2:m+1)' * Z0(1:n*n*n) - c;
      Z0(n*n*n+[1:2*m]) = [-min(z,0)+1e-3; max(z,0)+1e-3 ];

      % call SP
      [x,Z,ul,info,time] = sp(FF(:,1:m+1), blck_szs, c, ...
          x0, Z0, 10.0, abstol, reltol, 0.0, 100);

      % optimal value
      areas(i) = c'*x;

   end;

end;

figure(1);
ind = finite(areas);
plot(areas(ind), delays(ind));
xlabel('area');
ylabel('Tdom');



%%%%%%%%%%  two solutions  %%%%%%%%%%

% tolerances for SP
abstol = 1e-6;
reltol = 1e-6;

% points on the tradeoff curve
delays = [410, 2000];
sizes = zeros(m,2);

% compute tradeoff curve
for i=[1,2]

   delay = delays(i);

   disp([' ']);
   disp(['Compute solution ', int2str(i), ...
        ' (delay = ', num2str(delay), ').']);

   % n LMIs
   for j=1:n
      FF((j-1)*n*n+[1:n*n],1:m+1) = GG - (1/delay)*CC;
      FF((j-1)*n*n+(j-1)*n+j,1) = G0;
   end;
   
   % try very thick wire
   x0 = 0.99*wmax*ones(m,1); 
   t = Inf;
   for j=1:n
      t = min(t, min(eig(...
          reshape(FF((j-1)*n*n+[1:n*n],1:m+1)*[1;x0],n,n))));
   end;
   feasible = 1;

   if (t <= 0)   % x0 is infeasible, run phase 1

      disp([' ']);  disp(['Phase 1.']);

      % use slightly tighter bounds on x1, x2, d1, d2 to ensure
      % strict feasibility
      FF(n*n*n+[1:m],1) = -1e-5*ones(m,1);
      FF(n*n*n+m+[1:m],1) = (1-1e-5)*wmax*ones(m,1); 

      % dual feasible solution
      Z0 = zeros(n*n*n+2*m,1);
      for j=1:n
         Z0((j-1)*n*n+[1:n*n],1) = reshape(eye(n)/(n*n),n*n,1);
      end;
      z = FF(1:n*n*n,2:m+1)' * Z0(1:n*n*n);

      Z0(n*n*n+[1:2*m]) = [-min(z,0)+1e-3; max(z,0)+1e-3 ];

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
      FF(n*n*n+[1:m],1) = zeros(m,1);
      FF(n*n*n+m+[1:m],1) = wmax*ones(m,1); 

      % dual feasible point
      Z0 = zeros(n*n*n+2*m,1);
      for j=1:n
         Z0((j-1)*n*n+[1:n*n],1) = reshape(eye(n)/(n*n),n*n,1);
      end;
      z = FF(1:n*n*n,2:m+1)' * Z0(1:n*n*n) - c;
      Z0(n*n*n+[1:2*m]) = [-min(z,0)+1e-3; max(z,0)+1e-3 ];

      % call SP
      [sizes(:,i),Z,ul,info,time] = sp(FF(:,1:m+1), blck_szs, c, ...
          x0, Z0, 10.0, abstol, reltol, 0.0, 100);


   end;

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

 end;


keyboard

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

end;
