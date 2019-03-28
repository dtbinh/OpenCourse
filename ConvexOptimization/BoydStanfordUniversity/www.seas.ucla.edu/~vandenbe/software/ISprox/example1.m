%  Line spectra estimation by Toeplitz covariance fitting [CV18, Sec. IV-A]
%  This script calls the functions defined in ISproj.m and qmeif.m.
%  
%  A covariance matrix of the form
%                                        [     1       ][     1       ]'
%    R = sigma^2*I + sum_{k=1}^r |c_k|^2*[  e^(jw_k)   ][  e^(jw_k)   ]
%                                        [     :       ][     :       ]
%                                        [e^(j(n-1)w_k)][e^(j(n-1)w_k)]
%  is fit to a sample covariance matrix Rm by solving the semidefinite 
%  equivalence of the optimization problem 
%    minimize    rr*||R-Rm||_F^2 + sum_{k=1}^r |c_k|^2
%                                          [     1       ][     1       ]'
%    subject to  R = sigma^2*I             [  e^(jw_k)   ][  e^(jw_k)   ] 
%                    + sum_{k=1}^r |c_k|^2*[     :       ][     :       ]
%                                          [e^(j(n-1)w_k)][e^(j(n-1)w_k)]
%  with variables R, sigma^2, |c_k|, w_k, and r.
%  The parameters are then computed by doing a quadratic matrix 
%  factorization on the Toeplitz matrix R.
%  The sample covariance matrix is constructed as Rm = Y*Y'/m where Y is a 
%  n x m Hankel matrix with y(1),...,y(m) in its first row and y(m),...,
%  y(n+m-1) in its last column. The data series y of length N = n+m-1 
%  assume an underlying line spectrum model with circular white noise
%    y(i) = sum_{k=1}^r c_k*e^(ji*w_k) + v(i).
%  [CV16] H.-H. Chao and L. Vandenberghe , Semidefinite representations of 
%         gauge functions for structured low-rank matrix decomposition, 
%         2016. arXiv:1604.02500.
%  Author: Hsiao-Han (Suzi) Chao (suzihchao@ucla.edu)
%  Date modified: 7/18/2018

% change to use_octave = 1 if using Octave
use_octave = 0;

seed = 1;
rng(seed)
N = 150; % data length
amp = 5;
noise = 8; % standard deviation

%% Data generation
% ground truth
omegat =  pi/180 * [-18 -41 18 41]; % from [-pi,pi]
ct = amp * exp(1i*pi*rand(length(omegat),1)) .* randn(length(omegat),1); % coefficients
ct = [ct(1:2); conj(ct(1:2))];
yt = exp(1i * [0:N-1]' * omegat) * ct;
Rt = exp(1i * [0:N-1]' * omegat) * diag(ct.*conj(ct)) * exp(1i * [0:N-1]' * omegat)';

% adding noise
ym = yt + noise*randn(N,1); %(noise/sqrt(2)) * (randn(N,1) + 1i*randn(N,1));

figure(9);
plot(1:N,real(yt),'r','linewidth',2);
hold on;
plot(1:N,real(ym),'b:o','linewidth',1.5,'markersize',5);
hold off;

%% Recovery methods
n = 30;
m = 121;
rr = 0.0002;
% sampled covariance
Rm = (hankel(ym(1:n),ym(n:n+m-1))*hankel(ym(1:n),ym(n:n+m-1))')/m; 
DRm = zeros(n,1);
for k = 1:n
    DRm(k) = real(sum(diag(Rm,k-1)));
end
e0 = zeros(n,1); e0(1) = 1;
dD2 = (n:-1:1)';
dD = sqrt(dD2);

%% CVX
if ~use_octave
    cvx_begin sdp
        variable X(n,n) symmetric;
        variable x(n);
        minimize ((x-e0)'*(((1/(4*rr))*(x-e0)+DRm).*[1;2*ones(n-1,1)]./ (dD2)))
        subject to
            X >= 0;
            for k = 1:n
                x(k) == sum(diag(X,k-1));
            end
            x(1) == 1;
    cvx_end
else
    cvx_optval = -1.013119620642523e+002; % when seed = 1
end

%% Proximal gradient with IS distance 
p = n-1;
iMax = 2000; % maximum iteration
tol = 1e-4;
tolv = 1e-8; % tolerance for solving sub-problems
t = 10^(1)/p; % initial stepsize
beta = .5;
prox_time = cputime;
% initialization
xp = zeros(p+1,1); xp(1) = xp(1) + 1; % x = e_1
v = xp;
phiv = 0;
y = zeros(p+1,1);
gradphi_ = -xp;
theta = 1;
obj = ((xp-e0)'*(((1/(4*rr))*(xp-e0)+DRm).*[1;2*ones(n-1,1)]./ (dD2)));
obj_hist = zeros(iMax+1,1);
obj_hist(1) = obj;
gap_hist = zeros(iMax+1,1);
Stats = zeros(iMax,5); %[t; #iter; #newton; #levinson; cputime]
% iterate
for ii = 1:iMax
    theta = .5*(-theta^2+sqrt(theta^4+4*theta^2));
    %theta = 2/(ii+1); %
    %theta = 1; % no acceleration
    y = (1-theta)*xp + theta*v;
    grad_ = ((1/(2*rr))*(y-e0)+DRm)./dD2;
    u = -(t/theta)*grad_ + gradphi_;
    phivp = phiv;
    [vn,gradphi_n,phiv,stats] = ISproj(u,0,1,tolv);
        
    xp = (1-theta)*xp + theta*vn;
    objp = obj;
    obj = ((xp-e0)'*(((1/(4*rr))*(xp-e0)+DRm).*[1;2*ones(n-1,1)]./ (dD2)));
    % line search for t
    if t > 0
        rhs_tmp = (1-theta)*objp + ...
         theta*(((y-e0)'*(((1/(4*rr))*(y-e0)+DRm).*[1;2*ones(n-1,1)]./ (dD2))) ...
         - 2*(grad_(2:end)'*y(2:end))-grad_(1)*y(1));
        rhs = rhs_tmp + theta*(2*(grad_(2:end)'*vn(2:end))+grad_(1)*vn(1)+...
             (phiv-phivp-2*(gradphi_(2:end)'*(vn(2:end)-v(2:end)))...
             -gradphi_(1)*(vn(1)-v(1)))*theta/t);
        while obj > rhs*(1 + 2*p*eps)
            t = beta*t;
            disp([ii t]);
            disp(obj-rhs);
            u = -(t/theta)*grad_ + gradphi_;
            [vn,gradphi_n,phiv,stats] = ISproj(u,0,1,tolv);
            xp = (1-theta)*xp + theta*vn;
            obj = ((xp-e0)'*(((1/(4*rr))*(xp-e0)+DRm).*[1;2*ones(n-1,1)]./ (dD2)));
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

    % stopping condition
    if abs(obj-cvx_optval)/(-cvx_optval) < tol, break; end
end
prox_time = cputime - prox_time
if ii <= iMax
else
    ii = iMax;
end
Stats = Stats(1:ii,:);
obj_hist = obj_hist(1:ii+1);

figure; % convergence plot 
loglog(1:ii,(obj_hist(1:ii)-obj_hist(2:ii+1))./(-obj_hist(2:ii+1)),'m--','LineWidth',1.5);
hold on;
loglog(1:ii+1,(obj_hist(1:ii+1)-cvx_optval)/(-cvx_optval),'r','LineWidth',1.5);
xlabel('iteration','Fontsize', 20);
legend('(F(x^k) - F(x^{k-1}))/|F(x^{k-1})|','relative error against CVX');
set(gca, 'FontSize', 15);

% plot polynomial
figure;
ww = linspace(-pi,pi,6000)';
WW = xp(1) + 2*xp(2)*cos(ww);
for k = 2:p
    WW = WW + 2*xp(k+1)*cos(k*ww);
end
plot(ww,WW,'m-.','LineWidth',1.5);
hold on; 
plot([ww(1) ww(end)],[0 0],'k');
temp = axis; axis([-pi pi temp(3:4)]);

% recover supports and manitudes
dim = n;
Ty = toeplitz(((1/(2*rr))*(xp-e0)+DRm)./dD2);
[V, D] = eig(Ty);
tp = D(1,1);
d = diag(D);
d = d - tp*ones(n,1);
II  = find(d > 1e-4 * max(d));
r = length(II);
R = V(:,II) * diag(sqrt(d(II)));
[s,nu,~,~] = qmeif(R(1:dim-1,:),R(2:dim,:),[1 0;0 -1],0,1e-2);
omega = -angle(s./nu)';
Fc = exp(1i * [0:n-1]'*omega);
c = sqrt(diag(Fc \ (Fc\(Ty-tp*eye(n)))'));

% plot spectrum
stem(omega, abs(c),'b', 'LineWidth', 1, 'MarkerSize', 6); 
plot(omegat', abs(ct), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r')
xlabel('\omega', 'Fontsize', 20);
ylabel('magnitude', 'Fontsize', 20);
set(gca, 'FontSize', 15);
axis([-3, 3, 0, 1.2*amp]);
