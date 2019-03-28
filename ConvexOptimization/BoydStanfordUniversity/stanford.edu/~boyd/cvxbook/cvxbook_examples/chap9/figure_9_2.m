% Generates figure 9.2 of Boyd & Vandenberghe, Convex Optimization
%
% Gradient descent for minimizing f(x_1,x_2) = x_1^2 + 10*x_2^2.

M=10;

% plot some contour lines  
val=linspace(0,M*(M+1),5);  
theta=linspace(0,2*pi,100)';
plot(cos(theta)*sqrt(val), sin(theta)*sqrt(val/M), '--');
hold on;
axis('equal');

% 40 iterations
k=0:40;
plot(M*((M-1)/(M+1)).^k, (-(M-1)/(M+1)).^k,'-', ...
     M*((M-1)/(M+1)).^k, (-(M-1)/(M+1)).^k,'o');

set(gca,'XTick', [-10 0 10]);
set(gca,'YTick', [-4 0 4]);
xlabel('x1');
ylabel('x2');
axis([-12 12 -5 5]);
hold off
