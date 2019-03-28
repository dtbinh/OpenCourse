% Generates the linear dynamical system numerical example in the paper
% Sensor Selection via Convex Optimization
% www.stanford.edu/~boyd/papers/sensor_selection.html
%
% May 2008 Siddharth Joshi & Stephen Boyd
%
% Requires CVX; available from: http://www.stanford.edu/~boyd/cvx/

clear all;
randn('state',0); rand('state',0);
n = 5; p = 2; m = 100; q = m*p; k = 10;

% x(t+1) = Fx(t), t=1,...,m
% y(t) = Gx(t) + v(t), t=1,...,m
% dynamics matrix F
rho = [1.01 0.99]; 
theta = [0.03 0.05];
F = zeros(n,n);
F(1:2,1:2) = rho(1)*[cos(theta(1)) -sin(theta(1)); sin(theta(1)) cos(theta(1))];
F(3:4,3:4) = rho(2)*[cos(theta(2)) -sin(theta(2)); sin(theta(2)) cos(theta(2))];
F(5,5) = 0.98;

T = randn(n,n);
F = T*F*inv(T);
disp('Eigenvalues of F:'); disp(eig(F));

H = (2*rand(p,n) - 1)/10;
Sigma_1 = eye(n);
Fm= F^m;
Sigma_prior = Fm*Sigma_1*Fm';

A = zeros(q,n);
Finv = inv(F);
A(q-p+1:q,:) = H*Finv; 
for i=1:m-1
   A(q-p*(i+1)+1:q-p*i, :) = A(q-p*i+1:q-p*(i-1), :)*Finv;
end

cvx_begin 
	variable z(m)
    variable X(n,n) symmetric 
	minimize( trace(X) )
	subject to
        [X eye(n) ; eye(n) (A'*(kron(diag(z), eye(p)))*A + inv(Sigma_prior))] == semidefinite(2*n);
		z>=0;
		z<=1;
		sum(z) == k;
cvx_end

zsort=sort(z); thres=zsort(m-k); z01=(z>thres);

figure; hold on;
set(gca,'FontName','times', 'FontSize', 16);
xlabel('t'); ylabel('z');
plot(find(z01==1), 1, 'rx');
plot(find(z>0.001), z(z>0.001), 'bo'); 
hold off;

disp('Selected sensors:');
disp(find(z01==1))

mse01 = zeros(m,1); msecont = zeros(m,1);
I01 = inv(Sigma_prior);
Icont = I01;
for i = 1:m
    dyad = A(p*i-p+1:p*i, :)'*A(p*i-p+1:p*i, :);
    I01 = I01 + z01(i)*dyad;
    Icont = Icont + z(i)*dyad;
    mse01(i) = trace(inv(I01)); 
    msecont(i) = trace(inv(Icont));
end
figure; hold on;
set(gca,'FontName','times', 'FontSize', 16);
xlabel('t'); ylabel('mse');
stairs(0:m, [trace(Sigma_prior); mse01], 'b', 'LineWidth', 2); 
stairs(0:m,[trace(Sigma_prior); msecont], 'k:', 'LineWidth', 2);
hold off;
%print -deps ex2.eps
