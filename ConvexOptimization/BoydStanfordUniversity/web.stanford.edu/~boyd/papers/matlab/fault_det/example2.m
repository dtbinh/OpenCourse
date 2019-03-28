% EXAMPLE2.M
%
% Generates Figure 5 from
% "Relaxed Maximum a Posteriori Fault Identification"
% by A. Zymnis, S. Boyd and D. Gorinevsky
% requires CVX to run

clear all
cvx_quiet(false);

randn('state',32);
rand('state',34);


%%% ---------------- Generate Data -----------------------------------------
m = 10000; %number of sensors
n = 2000; %number of fault signatures
A = sprandn(m,n,0.01); %fault signatures
v1 = sum(A>0,1); v2 = sum(A>0,2);
A = A(setdiff(1:m,find(v2==0)),setdiff(1:n,find(v1==0)));
[m,n] = size(A);

pf = 0.002; %probability of fault
sigma = 1.5; %noise variance
lambda = log((1-pf)/pf)*ones(n,1);

x_true = double(rand(n,1)<pf); %true fault vector
y = A*x_true+sigma*randn(m,1); %measurement

%% --------------- Solve Relaxation ---------------------------------------
cvx_begin
    variable x(n)
    minimize((1/(2*sigma^2))*square_pos(norm(A*x-y,2))+lambda'*x)
    x >= 0
    x <= 1
cvx_end
l_min = cvx_optval-(1/(2*sigma^2))*norm(y)^2;

%% ------------- ARMAP with Local Opt. ----------------------------------
[X_amb_lc,l_amb_lc] = armap(A,y,pf*ones(n,1),sigma,0.01,1,'verbose');

%% -------------- Generate all Patterns by Enumeration --------------------
l_single = (1/(2*sigma^2))*norms(A-y*ones(1,n),2).^2+lambda'-(1/(2*sigma^2))*norm(y)^2*ones(1,n);

%% -------------- Generate Figures ----------------------------------------
figure; hold on; length = abs(l_min)/3;
line([l_min l_min],[0 length],'LineWidth',2,'Color','k');
line([0 0],[0 length],'LineWidth',2,'Color','k','LineStyle','--');
plot(l_amb_lc,length/3,'ro','MarkerSize',5);
plot(l_single(find(l_single<=0.2*abs(l_min))),2*length/3,'co','MarkerSize',5);
plot(l_amb_lc(1),length/3,'r.','MarkerSize',15);
plot(l_single(find(l_single==min(l_single))),2*length/3,'c.','MarkerSize',15);
xlabel('lyx')
axis equal
axis([l_min-0.1*abs(l_min) 0.2*abs(l_min) 0 length])
text(l_min,1.05*length,'llb')
text(0,1.05*length,'l0')
text(0.22*abs(l_min),length/3,'rmap')
text(0.22*abs(l_min),2*length/3,'single')
v = l_amb_lc(1);
text(v(1)-2,1*length/3,'xsrmap')
v = l_single(find(l_single==min(l_single)));
text(v(1)-2,2*length/3,'xsf')
set(gca,'YTick',[])
set(gca,'YColor',0.99*[1 1 1]);
set(gca,'FontSize',12)
set(gca,'Box','off')
%print -depsc penalty_ls.eps





