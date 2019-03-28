% CVX examples for ICM 2006 talk available at
% http://www.stanford.edu/~boyd/cvx_opt_graph_lapl_eigs.html
%
% This is an example with a larger graph (50 nodes and 200 edges).
% Written for CVX by Almir Mutapcic 08/29/06

%********************************************************************
% randomly generate a graph with 50 nodes and 200 edges
% and make it pretty for plotting
%********************************************************************
n = 50; threshold = 0.2529;
rand('state',209);
xy = rand(n,2);

angle = 10*pi/180;
Rotate = [ cos(angle) sin(angle); -sin(angle) cos(angle) ];
xy = (Rotate*xy')';

Dist = zeros(n,n);
for i=1:(n-1);
  for j=i+1:n;
    Dist(i,j) = norm( xy(i,:) - xy(j,:) );
  end;
end;
Dist = Dist + Dist';
Ad = Dist < threshold;
Ad = Ad - eye(n);
m = sum(sum(Ad))/2;

% find the incidence matrix
A = zeros(n,m);
l = 0;
for i=1:(n-1);
  for j=i+1:n;
    if Ad(i,j)>0.5
      l = l + 1;
      A(i,l) =  1;
      A(j,l) = -1;
    end;
  end;
end;
A = sparse(A);

%********************************************************************
% compute edge weights
%********************************************************************
fprintf(1,'WARNING: The optimal weight computations take some time...\n');
% size of the network
[n,m] = size(A);

% compute various edge weights (some optimal, some based on heuristics)
w_fdla = fdla(A);
w_fmmc = fmmc(A);
w_md   = max_deg(A);
w_bc   = best_const(A);
w_mh   = mh(A);

% report results
rho_fdla       = norm(eye(n) - A*diag(w_fdla)*A'  - 1/n*ones(n));
rho_fmmc       = norm(eye(n) - A*diag(w_fmmc)*A'  - 1/n*ones(n));
rho_max_deg    = norm(eye(n) - A*diag(w_md)*A'  - 1/n*ones(n));
rho_best_const = norm(eye(n) - A*diag(w_bc)*A'  - 1/n*ones(n));
rho_mh         = norm(eye(n) - A*diag(w_mh)*A'  - 1/n*ones(n));

tau_fdla       = 1/log(1/rho_fdla);
tau_fmmc       = 1/log(1/rho_fmmc);
tau_max_deg    = 1/log(1/rho_max_deg);
tau_best_const = 1/log(1/rho_best_const);
tau_mh         = 1/log(1/rho_mh);

eig_opt = sort(eig(eye(n) - A*diag(w_fdla)*A'));
eig_fmmc = sort(eig(eye(n) - A*diag(w_fmmc)*A'));
eig_mh = sort(eig(eye(n) - A*diag(w_mh)*A'));
eig_md = sort(eig(eye(n) - A*diag(w_md)*A'));
eig_bc = sort(eig(eye(n) - A*diag(w_bc)*A'));

fprintf(1,'\nResults:\n');
fprintf(1,'FDLA weights:\t\t rho = %5.4f \t tau = %5.4f\n',rho_fdla,tau_fdla);
fprintf(1,'FMMC weights:\t\t rho = %5.4f \t tau = %5.4f\n',rho_fmmc,tau_fmmc);
fprintf(1,'M-H weights:\t\t rho = %5.4f \t tau = %5.4f\n',rho_mh,tau_mh);
fprintf(1,'MAX_DEG weights:\t rho = %5.4f \t tau = %5.4f\n',rho_max_deg,tau_max_deg);
fprintf(1,'BEST_CONST weights:\t rho = %5.4f \t tau = %5.4f\n', ...
        rho_best_const,tau_best_const);

%********************************************************************
% plot results
%********************************************************************
figure(1)
gplot(Ad,xy);
hold on;
plot(xy(:,1), xy(:,2), 'ko','LineWidth',4, 'MarkerSize',4);
axis([0.05 1.1 -0.1 0.95]);
title('Graph')
hold off;

figure(2)
v_fdla = [w_fdla; diag(eye(n) - A*diag(w_fdla)*A')];
[ifdla, jfdla, neg_fdla] = find( v_fdla.*(v_fdla < -0.001 ) );
v_fdla(ifdla) = [];
wbins = [-0.6:0.012:0.6];
hist(neg_fdla,wbins); hold on,
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r')
hist(v_fdla,wbins); hold off,
axis([-0.6 0.6 0 12]);
xlabel('optimal FDLA weights');
ylabel('histogram');

figure(3)
xbins = (-1:0.015:1)';
ymax  = 6;
subplot(3,1,1)
hist(eig_md, xbins); hold on;
max_md = max(abs(eig_md(1:n-1)));
plot([-max_md -max_md],[0 ymax], 'b--');
plot([ max_md  max_md],[0 ymax], 'b--');
axis([-1 1 0 ymax]);
text(0,5,'MAX DEG');
title('Eigenvalue distributions')
subplot(3,1,2)
hist(eig_bc, xbins); hold on;
max_opt = max(abs(eig_bc(1:n-1)));
plot([-max_opt -max_opt],[0 ymax], 'b--');
plot([ max_opt  max_opt],[0 ymax], 'b--');
axis([-1 1 0 ymax]);
text(0,5,'BEST CONST');
subplot(3,1,3)
hist(eig_opt, xbins); hold on;
max_opt = max(abs(eig_opt(1:n-1)));
plot([-max_opt -max_opt],[0 ymax], 'b--');
plot([ max_opt  max_opt],[0 ymax], 'b--');
axis([-1 1 0 ymax]);
text(0,5,'FDLA');

figure(4)
xbins = (-1:0.015:1)';
ymax  = 6;
subplot(3,1,1)
hist(eig_md, xbins); hold on;
max_md = max(abs(eig_md(1:n-1)));
plot([-max_md -max_md],[0 ymax], 'b--');
plot([ max_md  max_md],[0 ymax], 'b--');
axis([-1 1 0 ymax]);
text(0,5,'MAX DEG');
title('Eigenvalue distributions')
subplot(3,1,2)
hist(eig_mh, xbins); hold on;
max_opt = max(abs(eig_mh(1:n-1)));
plot([-max_opt -max_opt],[0 ymax], 'b--');
plot([ max_opt  max_opt],[0 ymax], 'b--');
axis([-1 1 0 ymax]);
text(0,5,'MH');
subplot(3,1,3)
hist(eig_fmmc, xbins); hold on;
max_opt = max(abs(eig_fmmc(1:n-1)));
plot([-max_opt -max_opt],[0 ymax], 'b--');
plot([ max_opt  max_opt],[0 ymax], 'b--');
axis([-1 1 0 ymax]);
text(0,5,'FMMC');

figure(5)
v_fmmc = [w_fmmc; diag(eye(n) - A*diag(w_fmmc)*A')];
[ifmmc, jfmmc, nonzero_fmmc] = find( v_fmmc.*(v_fmmc > 0.001 ) );
hist(nonzero_fmmc,80);
axis([0 1 0 10]);
xlabel('optimal positive FMMC weights');
ylabel('histogram');

figure(6)
An = abs(A*diag(w_fmmc)*A');
An = (An - diag(diag(An))) > 0.0001;
gplot(An,xy,'b-'); hold on;
h = findobj(gca,'Type','line');
set(h,'LineWidth',2.5)
gplot(Ad,xy,'b:');
plot(xy(:,1), xy(:,2), 'ko','LineWidth',4, 'MarkerSize',4);
axis([0.05 1.1 -0.1 0.95]);
title('Subgraph with positive transition prob.')
hold off;
