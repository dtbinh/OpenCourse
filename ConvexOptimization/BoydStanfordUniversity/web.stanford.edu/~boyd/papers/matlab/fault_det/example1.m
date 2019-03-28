% EXAMPLE1.M
%
% Generates Figure 1 from
% "Relaxed Maximum a Posteriori Fault Identification"
% by A. Zymnis, S. Boyd and D. Gorinevsky
% requires CVX to run

clear all
cvx_quiet(true);

randn('state',824);
rand('state',25);


%% ---------------- Generate Data -----------------------------------------
m = 50; %number of sensors
n = 100; %number of fault signatures
pf = 0.05; %probability of fault
sigma = 1; %noise std
lambda = log((1-pf)/pf)*ones(n,1);
K = 10; %ambiguity set size


A = randn(m,n); %fault signatures
x_true = double(rand(n,1)<pf); %true fault vector
X_true = x_true*ones(1,K);
y = A*x_true+sigma*randn(m,1); %measurement
l_true = (1/(2*sigma^2))*norm(A*x_true-y)^2+lambda'*x_true-(1/(2*sigma^2))*norm(y)^2;

%% -------------- Get Lower Bound Using CVX -------------------------------------
cvx_begin
    variable x(n)
    minimize((1/(2*sigma^2))*square_pos(norm(A*x-y,2))+lambda'*x)
    x >= 0
    x <= 1
cvx_end
l_min = cvx_optval-(1/(2*sigma^2))*norm(y)^2;

%% -------------- ARMAP Without Local Opt ------------------------------------
[X_amb,l_amb] = armap(A,y,pf*ones(n,1),sigma,1/(2*n),0,'verbose');

%% -------------- ARMAP With Local Opt ---------------------------------------
[X_amb_lc,l_amb_lc] = armap(A,y,pf*ones(n,1),sigma,1/(2*n),1,'verbose');

%% -------------- Generate all Patterns by Enumeration --------------------
fprintf(1,'Generating 1-fault and 2-fault patterns...\n')
x_single = eye(n); 
l_single = (1/(2*sigma^2))*norms(A*x_single-y*ones(1,n),2).^2+lambda'*x_single-(1/(2*sigma^2))*norm(y)^2*ones(1,n);
[l_single,ind] = sort(l_single,'ascend');

x_double = []; l_double = [];
for i = 1:n-1
    for j = i+1:n
        x_cur = zeros(n,1);
        x_cur(i) = 1; x_cur(j) = 1;
        l_cur = (1/(2*sigma^2))*norm(A*x_cur-y)^2+lambda'*x_cur-(1/(2*sigma^2))*norm(y)^2;
        x_double = [x_double x_cur]; l_double = [l_double l_cur];
    end
end
[l_double,ind] = sort(l_double,'ascend');

%% -------------- Generate Figure -----------------------------------------
figure; hold on; 
v_max = 0.2*abs(l_min);
length = (v_max+abs(l_min))/2;
line([l_min l_min],[0 length],'LineWidth',2,'Color','k');
line([0 0],[0 length],'LineWidth',2,'Color','k','LineStyle','--');
plot(l_amb(1:K),2*length/5,'ro','MarkerSize',5);
plot(l_amb_lc,1*length/5,'mo','MarkerSize',5);
plot(l_single(1:K),4*length/5,'co','MarkerSize',5);
plot(l_double(1:K),3*length/5,'bo','MarkerSize',5);
plot(l_amb(1),2*length/5,'r.','MarkerSize',15);
plot(l_amb_lc(1),1*length/5,'m.','MarkerSize',15);
plot(l_single(1),4*length/5,'c.','MarkerSize',15);
plot(l_double(1),3*length/5,'b.','MarkerSize',15);
xlabel('lyx')
axis equal
axis([l_min-0.1*abs(l_min) v_max 0 length])
text(l_min,1.05*length,'llb')
text(0,1.05*length,'l0')
text(1.05*v_max,2*length/5,'rmap')
text(1.05*v_max,1*length/5,'locopt')
text(1.05*v_max,4*length/5,'single')
text(1.05*v_max,3*length/5,'double')
v = l_amb(1);
text(v(1)-0.5,2*length/5,'xsrmap')
v = l_amb_lc(1);
text(v(1)-0.5,1*length/5,'xrrmap')
v = l_single(find(l_single==min(l_single)));
text(v(1)-0.5,4*length/5,'xsf')
v = l_double(find(l_double==min(l_double)));
text(v(1)-0.5,3*length/5,'xdf')
set(gca,'YTick',[])
set(gca,'YColor',0.99*[1 1 1]);
set(gca,'FontSize',12)
set(gca,'Box','off')
print -depsc penalty.eps

return
