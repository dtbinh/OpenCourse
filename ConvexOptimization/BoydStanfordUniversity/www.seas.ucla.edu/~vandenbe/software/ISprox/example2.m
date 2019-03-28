%  This script tests accelerated proximal gradient algorithm with Bregman 
%  distance and Itakura-Saito kernel function on the optimization problem
%  [CV18, Sec. IV-B]
%  min  ||x-b||^2
%  s.t. X(w) = x_0 + 2 sum_{k=1}^p x_k cos(kw) >= 0,  for 0 <= w < 2\pi
%      x_0 = 1.
%
%  This script calls the function defined in ISproj.m.
%
%  Author: Hsiao-Han (Suzi) Chao (suzihchao@ucla.edu)
%  Date modified: 7/18/2018

if ~exist('benching','var')
use_octave = 0; % change to use_octave = 1 when using Octave
use_plot = 1;

% problem setup
seed = 509;
rng(seed);
p = 399;
b = randn(p+1,1); 
end

if use_plot
lineW = 2;
% plot b_0 + 2 sum_{k=1}^p b_k cos(kw)
ww = linspace(0,2*pi,1000)';
WW = b(1) + 2*b(2)*cos(ww);
for k = 2:p
    WW = WW + 2*b(k+1)*cos(k*ww);
end
figure; plot(ww,WW,'r:','LineWidth',lineW); hold on;
end

if ~use_octave
%% SDP formulation solved via CVX
% min  ||x-b||^2
% s.t. x_k = tr E_k^T Z, k = 0,1,...,p
%      Z >= 0

fprintf('   CVX is being called...\n');
%cvx_solver sedumi
cvx_begin sdp 
    %cvx_precision('high')
    variable x(p+1);
    variable Z(p+1,p+1) symmetric;
    minimize ((x-b)'*(x-b));
    subject to
        for k = 0:p
            x(k+1) == sum(diag(Z,k));
        end
        Z >= 0;
        x(1) == 1;
cvx_end
    
x_cvx = x;

fprintf('   cputime for CVX (seconds): %f \n\n',cvx_cputime);

if use_plot
% plot X(w)
WW = x(1) + 2*x(2)*cos(ww);
for k = 2:p
    WW = WW + 2*x(k+1)*cos(k*ww);
end
plot(ww,WW,'b','LineWidth',lineW);
temp = axis; axis([0 2*pi temp(3:4)]);
end
end

%% Proximal gradient with IS distance
fprintf('   The projected gradient method is started...\n');
iMax = 1000; % maximum iteration
if use_octave
    tol = 1e-6;
else
    tol = 1e-4;
end
tolv = 1e-8; % tolerance for solving sub-problems
t = 10^(1)/p; % initial stepsize
beta = .5;
prox_time = cputime;
% initialization
x = zeros(p+1,1); x(1) = x(1) + 1; % x = e_1
v = x;
phiv = 0;
y = zeros(p+1,1);
gradphi_ = -x;
theta = 1;
obj = norm(x-b)^2;
objp = obj;
obj_hist = zeros(iMax+1,1);
obj_hist(1) = obj;
if ~use_octave
    err_hist = zeros(iMax+1,1);
    err_hist(1) = norm(x-x_cvx);
end
Stats = zeros(iMax,5); %[t; #iter; #newton; #levinson; cputime]
% iterate
for ii = 1:iMax
    theta = .5*(-theta^2+sqrt(theta^4+4*theta^2));
    %theta = 2/(ii+1); %
    %theta = 1; % no acceleration
    y = (1-theta)*x + theta*v;
    grad_ = y-b; % 1/2 gradient
    grad_(1) = grad_(1)*2;
    u = -(t/theta)*grad_ + gradphi_;
    phivp = phiv;
    [vn,gradphi_n,phiv,stats] = ISproj(u,0,1,tolv);
        
    x = (1-theta)*x + theta*vn;
    objp = min(obj,objp);
    obj = norm(x-b)^2;
    % line search for t
    if t > 0
    rhs_tmp = (1-theta)*objp + ...
     theta*(norm(y-b)^2 - 2*(grad_(2:end)'*y(2:end))-grad_(1)*y(1));
    rhs = rhs_tmp + theta*(2*(grad_(2:end)'*vn(2:end))+grad_(1)*vn(1)+...
         (phiv-phivp-2*(gradphi_(2:end)'*(vn(2:end)-v(2:end)))...
         -gradphi_(1)*(vn(1)-v(1)))*theta/t);
  
    while obj > rhs*(1 + 2*p*eps)
     t = beta*t;
     u = -(t/theta)*grad_ + gradphi_;
     [vn,gradphi_n,phiv,stats] = ISproj(u,0,1,tolv);
     x = (1-theta)*x + theta*vn;
     obj = norm(x-b)^2;
     rhs = rhs_tmp + theta*(2*(grad_(2:end)'*vn(2:end))+grad_(1)*vn(1)+...
         (phiv-phivp-2*(gradphi_(2:end)'*(vn(2:end)-v(2:end)))...
         -gradphi_(1)*(vn(1)-v(1)))*theta/t);
    end
    end
    v = vn;
    gradphi_n(2:end) = gradphi_n(2:end)/2;
    gradphi_ = gradphi_n;

    obj_hist(ii+1) = obj;
    Stats(ii,:) = [t stats'];
    if ~use_octave
        err_hist(ii+1) = norm(x-x_cvx);
    end

    % stopping condition
    if use_octave
        if abs(obj-objp)/min(obj,objp) < tol, break; end
    else
        if abs(obj-cvx_optval)/cvx_optval < tol, break; end
    end
end
prox_time = cputime - prox_time;
fprintf('   cputime for projected gradient method (seconds): %f \n',prox_time);
if ii <= iMax
else
    ii = iMax;
end
Stats = Stats(1:ii,:);
obj_hist = obj_hist(1:ii+1);

if use_plot
    % plot X(w)
    WW = x(1) + 2*x(2)*cos(ww);
    for k = 2:p
        WW = WW + 2*x(k+1)*cos(k*ww);
    end
    plot(ww,WW,'m-.','LineWidth',lineW);
    plot([ww(1) ww(end)],[0 0],'k');
    temp = axis; axis([0 2*pi temp(3:4)]);
    if use_octave
        legend('B(w)','X(w) from prox');
    else
        legend('F_a(e^{jw})','F_x(e^{jw}) from CVX','F_x(e^{jw}) from prox');
    end

    if use_octave
        figure; loglog(1:ii,obj_hist(1:ii)-obj_hist(2:ii+1));
        ylabel('F(x^k) - F(x^{k-1})');
    else
        figure; loglog(1:ii+1,obj_hist(1:ii+1)-cvx_optval);
        hold on; loglog(1:ii+1,err_hist(1:ii+1),'r');
        legend('F(x^k) - cvx\_optval','||x-x\_cvx||')
    end
    xlabel('iteration');
    title(['p = ' num2str(p) ', rng(' num2str(seed) ')']);
end