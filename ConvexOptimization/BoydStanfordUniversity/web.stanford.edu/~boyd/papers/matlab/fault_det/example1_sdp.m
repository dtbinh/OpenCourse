% EXAMPLE1_SDP.M
%
% Compares QP and SDP based relaxation methods for problem
% described in paper:
% "Relaxed Maximum a Posteriori Fault Identification"
% by A. Zymnis, S. Boyd and D. Gorinevsky
% requires CVX to run
% WARNING: Augmented SDP takes a while to run

clear all
cvx_quiet(false);

randn('state',824);
rand('state',25);


%% ---------------- Generate Data -----------------------------------------
m = 50; %number of sensors
n = 100; %number of fault signatures
pf = 0.05; %probability of fault
sigma = 1.5; %noise std
lambda = log((1-pf)/pf)*ones(n,1);
K = 10; %ambiguity set size


A = randn(m,n); %fault signatures
x_true = double(rand(n,1)<pf); %true fault vector
X_true = x_true*ones(1,K);
y = A*x_true+sigma*randn(m,1); %measurement
l_true = (1/(2*sigma^2))*norm(A*x_true-y)^2+lambda'*x_true-(1/(2*sigma^2))*norm(y)^2;

%% -------------- Solve using QP relaxation -------------------------------------

% solve relaxation
cvx_begin
    variable x(n)
    minimize((1/(2*sigma^2))*square_pos(norm(A*x-y,2))+lambda'*x)
    x >= 0
    x <= 1
cvx_end
l_min_qp = cvx_optval-(1/(2*sigma^2))*norm(y)^2

%% rounding
[x_sort,ind_x] = sort(x,'descend'); x_cand = []; l_cand = [];
for i = 1:n
    x_cur = zeros(n,1);
    x_cur(ind_x(1:i)) = 1;
    l_cur = (1/(2*sigma^2))*square_pos(norm(A*x_cur-y,2))+lambda'*x_cur-(1/(2*sigma^2))*norm(y)^2;
    x_cand = [x_cand x_cur]; l_cand = [l_cand l_cur]; 
end
[l_sort,ind_l] = sort(l_cand,'ascend');
X_amb_qp = x_cand(:,ind_l(1:K)); %get ambiguity set
l_amb_qp = l_sort(1:K);

% perform local optimization
EXIT_FLAG = 0; iter = 0;
while(~EXIT_FLAG)
    x_cur = X_amb_qp(:,1); x_best = x_cur;
    for i = 1:n
        iter = iter+1;
        x_cur(i) = not(x_cur(i));
        l_cur = (1/(2*sigma^2))*norm(A*x_cur-y,2).^2+lambda'*x_cur-(1/(2*sigma^2))*norm(y)^2;
        if any(l_cur<l_amb_qp)
            ind = find(l_cur<l_amb_qp);
            ind = ind(1);
            X_amb_qp = [X_amb_qp(:,1:(ind-1)) x_cur X_amb_qp(:,ind:(end-1))];
            l_amb_qp = [l_amb_qp(:,1:(ind-1)) l_cur l_amb_qp(:,ind:(end-1))];
            if ind==1
                fprintf(1,'Found new best pattern!\n');
            else
                x_cur(i) = not(x_cur(i));
            end
        else
            x_cur(i) = not(x_cur(i));
        end
    end
    if all(x_best == X_amb_qp(:,1)), 
        EXIT_FLAG = 1; 
    end
end

l_amb_qp


%% -------------- Solve using SDP relaxation ------------------

% solve SDP relaxation
W = (1/(2*sigma^2))*A'*A;
w = lambda-(1/sigma^2)*A'*y;
mu = (1/(2*sigma^2))*y'*y;
cvx_begin
    variable z(n)
    variable Z(n,n) symmetric
    minimize(trace(W*Z)+w'*z+mu)
    subject to
        diag(Z)==z;
        [Z z; z' 1] == semidefinite(n+1);
        z >= 0;
        z <= 1;
cvx_end
l_min_sdp = cvx_optval-(1/(2*sigma^2))*norm(y)^2

%% rounding
[x_sort,ind_x] = sort(z,'descend'); x_cand = []; l_cand = [];
for i = 1:n
    x_cur = zeros(n,1);
    x_cur(ind_x(1:i)) = 1;
    l_cur = (1/(2*sigma^2))*square_pos(norm(A*x_cur-y,2))+lambda'*x_cur-(1/(2*sigma^2))*norm(y)^2;
    x_cand = [x_cand x_cur]; l_cand = [l_cand l_cur]; 
end
[l_sort,ind_l] = sort(l_cand,'ascend');
X_amb_sdp = x_cand(:,ind_l(1:K)); %get ambiguity set
l_amb_sdp = l_sort(1:K);


% perform local optimization
EXIT_FLAG = 0; iter = 0;
while(~EXIT_FLAG)
    x_cur = X_amb_sdp(:,1); x_best = x_cur;
    for i = 1:n
        iter = iter+1;
        x_cur(i) = not(x_cur(i));
        l_cur = (1/(2*sigma^2))*norm(A*x_cur-y,2).^2+lambda'*x_cur-(1/(2*sigma^2))*norm(y)^2;
        if any(l_cur<l_amb_sdp)
            ind = find(l_cur<l_amb_sdp);
            ind = ind(1);
            X_amb_sdp = [X_amb_sdp(:,1:(ind-1)) x_cur X_amb_sdp(:,ind:(end-1))];
            l_amb_sdp = [l_amb_sdp(:,1:(ind-1)) l_cur l_amb_sdp(:,ind:(end-1))];
            if ind==1
                fprintf(1,'Found new best pattern!\n');
            else
                x_cur(i) = not(x_cur(i));
            end
        else
            x_cur(i) = not(x_cur(i));
        end
    end
    if all(x_best == X_amb_sdp(:,1)), 
        EXIT_FLAG = 1; 
    end
end

l_amb_sdp1 = l_amb_sdp


% solve tighter SDP relaxation
W = (1/(2*sigma^2))*A'*A;
w = lambda-(1/sigma^2)*A'*y;
mu = (1/(2*sigma^2))*y'*y;
cvx_begin
    variable z(n)
    variable Z(n,n) symmetric
    minimize(trace(W*Z)+w'*z+mu)
    subject to
        diag(Z)==z;
        [Z z; z' 1] == semidefinite(n+1);
        Z >= 0;
        z >= 0;
        z <= 1;
cvx_end
l_min_sdp2 = cvx_optval-(1/(2*sigma^2))*norm(y)^2

%% rounding
[x_sort,ind_x] = sort(z,'descend'); x_cand = []; l_cand = [];
for i = 1:n
    x_cur = zeros(n,1);
    x_cur(ind_x(1:i)) = 1;
    l_cur = (1/(2*sigma^2))*square_pos(norm(A*x_cur-y,2))+lambda'*x_cur-(1/(2*sigma^2))*norm(y)^2;
    x_cand = [x_cand x_cur]; l_cand = [l_cand l_cur]; 
end
[l_sort,ind_l] = sort(l_cand,'ascend');
X_amb_sdp = x_cand(:,ind_l(1:K)); %get ambiguity set
l_amb_sdp = l_sort(1:K);

% perform local optimization
EXIT_FLAG = 0; iter = 0;
while(~EXIT_FLAG)
    x_cur = X_amb_sdp(:,1); x_best = x_cur;
    for i = 1:n
        iter = iter+1;
        x_cur(i) = not(x_cur(i));
        l_cur = (1/(2*sigma^2))*norm(A*x_cur-y,2).^2+lambda'*x_cur-(1/(2*sigma^2))*norm(y)^2;
        if any(l_cur<l_amb_sdp)
            ind = find(l_cur<l_amb_sdp);
            ind = ind(1);
            X_amb_sdp = [X_amb_sdp(:,1:(ind-1)) x_cur X_amb_sdp(:,ind:(end-1))];
            l_amb_sdp = [l_amb_sdp(:,1:(ind-1)) l_cur l_amb_sdp(:,ind:(end-1))];
            if ind==1
                fprintf(1,'Found new best pattern!\n');
            else
                x_cur(i) = not(x_cur(i));
            end
        else
            x_cur(i) = not(x_cur(i));
        end
    end
    if all(x_best == X_amb_sdp(:,1)), 
        EXIT_FLAG = 1; 
    end
end

l_amb_sdp2 = l_amb_sdp


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
v_max = 0.01*abs(l_min_qp);
length = (v_max+abs(l_min_qp))/2;
line([l_min_qp l_min_qp],[0 length],'LineWidth',2,'Color','k');
line([l_min_sdp l_min_sdp],[0 length],'LineWidth',2,'Color','b');
line([l_min_sdp2 l_min_sdp2],[0 length],'LineWidth',2,'Color','r');
line([0 0],[0 length],'LineWidth',2,'Color','k','LineStyle','--');
plot(l_amb_qp,3*length/6,'ko','MarkerSize',5);
plot(l_amb_sdp2(1:K),1*length/6,'ro','MarkerSize',5);
plot(l_amb_sdp1(1:K),2*length/6,'bo','MarkerSize',5);
plot(l_single(1:K),5*length/6,'mo','MarkerSize',5);
plot(l_double(1:K),4*length/6,'co','MarkerSize',5);
plot(l_amb_qp(1),3*length/6,'k.','MarkerSize',15);
plot(l_amb_sdp2(1),1*length/6,'r.','MarkerSize',15);
plot(l_amb_sdp1(1),2*length/6,'b.','MarkerSize',15);
plot(l_single(1),5*length/6,'m.','MarkerSize',15);
plot(l_double(1),4*length/6,'c.','MarkerSize',15);
xlabel('lyx')
axis equal
axis([l_min_qp-0.05*abs(l_min_qp) v_max 0 length])
text(l_min_qp,1.05*length,'llb')
text(l_min_sdp,1.05*length,'lsdp')
text(l_min_sdp2,1.05*length,'lsdp2')
text(0,1.05*length,'l0')
text(1.05*v_max,3*length/6,'rmap')
text(1.05*v_max,1*length/6,'sdp2')
text(1.05*v_max,2*length/6,'sdp1')
text(1.05*v_max,5*length/6,'single')
text(1.05*v_max,4*length/6,'double')
v = l_amb_qp(1);
text(v(1)+0.5,1*length/6,'xrmap')
v = l_amb_sdp1(1);
text(v(1)+0.5,2*length/6,'xsdp1')
v = l_amb_sdp2(1);
text(v(1)+0.5,3*length/6,'xsdp2')
v = l_single(find(l_single==min(l_single)));
text(v(1)+0.5,5*length/6,'xsf')
v = l_double(find(l_double==min(l_double)));
text(v(1)+0.5,4*length/6,'xdf')
set(gca,'YTick',[])
set(gca,'YColor',0.99*[1 1 1]);
set(gca,'FontSize',12)
set(gca,'Box','off')
print -depsc penalty_sdp2.eps



