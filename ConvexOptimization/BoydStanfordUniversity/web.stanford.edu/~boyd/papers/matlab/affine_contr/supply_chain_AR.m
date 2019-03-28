function [F, u0, x_tilde, u_tilde, Pxw, Puw, cvx_status] = supply_chain_AR(x_0, c, b, w, wf, mu, Sigma, M)
% SUPPLY_CHAIN_AR computes the affine controller that minimizes
%		  sum_{t=0}^T(Phi(x(t)) subject to 
% 	  	  x(t+1) = x(t) + u(t) - d(t) 
% where `d ~ LN(mu, Sigma) (here d is the whole trajectory).
% Here T is the length of the trajectory d (i.e., the length 
% of mu). 

T = length(mu); 

G = ones(T+1);
G = tril(G) - diag(diag(G)); 
G = G(:,1:end-1); 
H = G; 
x0 = x_0*ones(T+1,1); 

% optimal affine controller: 
W = -exp(mu*ones(1,M) + sqrtm(Sigma)*randn(T,M));
I = eye(T+1); 
cvx_begin
    variable Q(T,T+1) 
    variable r(T) 
       
    Pxw = (I+H*Q)*G;
    Puw = Q*G;
    x_tilde = (I + H*Q)*x0 + H*r;
    u_tilde = Q*x0 + r;
    
    x = Pxw*W + x_tilde*ones(1,M);
    u = Puw*W + u_tilde*ones(1,M);
    minimize(sum(sum(max(w*x(1:end-1,:),-b*x(1:end-1,:))) + max(wf*x(end,:),-b*x(end,:)) + sum(abs(c*u))))
    Qupper = [Q;zeros(1,T+1)];
    Qupper = triu(Qupper) - diag(diag(Qupper));
    Qupper == 0
cvx_end

I = eye(T);
F = (I+Q*H)\Q;
u0 = (I+Q*H)\r;
