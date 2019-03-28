function [X_amb, l_amb, ERROR_FLAG] = armap(A,y,p,sigma,kappa,LO_FLAG,verbosity)
%
% [X_amb, l_amb, ERROR_FLAG] = ARMAP(A,y,p,sigma,kappa,LO_FLAG,verbosity)
%
% Solves the approximate relaxed MAP fault identification problem:
% 
% minimize l_y(x) + kappa*psi(x)
%
% This problem is described in
% "Relaxed Maximum a Posteriori Fault Identification"
% by A. Zymnis, S. Boyd, and D. Gorinevsky
%
% Input parameters:
% A: m by n matrix of fault signatures
% y: m by 1 vector of continuous measurements 
%    OR m by 2 matrix of quantized intervals
% p: n vector of prior fault probabilities
% sigma: noise standard deviation
% kappa: log-barrier penalty parameter
% LO_FLAG: (optional) perform local optimization if set to 1
%          default is 1
% verbosity: (optional) if equal to 'quiet', no text output
%            default is verbose
%
% Returned variables:
% X_amb: n by K matrix whose kth column is the 
%        kth element in the ambiguity group
% l_amb: K vector of losses in ambiguity group
% ERROR_FLAG: 0, if everything OK
%             1, if numerical errors occured
%

%% Problem setup

% Problem data
ERROR_FLAG = 0; X_amb = []; l_amb = [];
[m,n] = size(A);
lambda = log((1-p)./p);
QUANT_FLAG = 0;
if size(y,2) == 2
    t_min = y(:,1);
    t_max = y(:,2);
    y = (t_min+t_max)/2;
    QUANT_FLAG = 1;
end
warning off

% Newton method parameters
alpha = 0.0001;
beta = 0.5;
MAX_ITER = 50;
tol = 1e-5;

% Local optimization parameters 
K = 10; %size of ambiguity set
num_sort = max(K,3*sum(p)); %number of sorted fault patterns
LO_ITER = 10*n; %maximum number of local optimization iterations

% Functions that generate gradient and Hessian
ATA = A'*A; ATy = A'*y;
g_pen = @(x)(1./x-1./(1-x));
H_pen = @(x)(spdiags(1./x.^2,0,n,n)+...
    spdiags(1./(1-x).^2,0,n,n));
g_cont = @(x)((1/sigma^2)*(ATA*x - ATy)+lambda);
H_cont = @(x)((1/sigma^2)*ATA);
g_quant = @(l,u)(lambda-(1/sigma)*A'*((exp(-l.^2/2)-exp(-u.^2/2))./...
    (sqrt(2*pi)*(normcdf(u)-normcdf(l)))));
H_quant = @(l,u)((1/sigma^2)*A'*spdiags(...
    (-l.*exp(-l.^2/2)+u.*exp(-u.^2/2))./(sqrt(2*pi)*(normcdf(u)-normcdf(l)))+...
    ((exp(-l.^2/2)-exp(-u.^2/2))./(sqrt(2*pi)*(normcdf(u)-normcdf(l)))).^2,0,m,m)*A);

% Function that generates fault pattern loss
if ~QUANT_FLAG
    l_0 = (1/(2*sigma^2))*norm(y).^2;
    l_y = @(x)((1/(2*sigma^2))*norm(A*x-y).^2+lambda'*x-l_0);
else
    l_0 = -sum(log(normcdf(t_max/sigma)-normcdf(t_min/sigma)));
    l_y = @(x)(-sum(log(normcdf((-A*x+t_max)/sigma)-...
                    normcdf((-A*x+t_min)/sigma)))+lambda'*x -l_0);
end

% Set verbosity
if nargin<=5
    LO_FLAG = 1;
end
if nargin<=6
    verb = 1;
else
    verb = ~strcmp(verbosity,'quiet');
end

if verb
    fprintf(1,'\nFault detection problem with n = %d possible faults and m = %d sensors.\n\n',...
        n,m);
end

%% Newton method for continuous measurements
x = 0.5*ones(n,1);
if verb
    if ~QUANT_FLAG 
        fprintf(1,'Performing Newton method...\n');
    else 
        fprintf(1,'Finding feasible point...\n'); 
    end
end
for i = 1:MAX_ITER
    if i == MAX_ITER, ERROR_FLAG = 1; end
    
    % Form gradient and Hessian
    g = g_cont(x)-kappa*g_pen(x);
    H = H_cont(x)+kappa*H_pen(x);
    
    % Compute Newton step
    dx = H\(-g);
    
    % Get into feasible region
    gamma = 1;
    if max(x+gamma*dx)>1
        ind_x = find(dx>0);
        gamma = min(gamma,0.99*min((1-x(ind_x))./dx(ind_x)));
    end
    
    if min(x+gamma*dx)<0
        ind_x = find(dx<0);
        gamma = min(gamma,0.99*min(-x(ind_x)./dx(ind_x)));
    end
    
    % Line search
    norm_cur = norm(g);
    while (1)
        g_next = g_cont(x+gamma*dx)-kappa*g_pen(x+gamma*dx);
        norm_next = norm(g_next);
        if norm_next <= (1-alpha*gamma)*norm_cur, break; end
        gamma = beta*gamma;
    end
    
    x=x+gamma*dx;
    
    norm_ratio = norm_next/sqrt(n);
    if verb
        fprintf(1,'Iteration: %2d, Step size: %e, Loss: %e\n',i,gamma,l_y(x))
    end
    
    if norm_ratio<1e-3, break; end
end
if verb, fprintf(1,'Done!\n\n'); end

if ERROR_FLAG==1
    if verb, fprintf(1,'Ran into numerical problems...\n\n'); end
    return
end
        

%% Newton method for quantized measurements
if QUANT_FLAG
if verb
    fprintf(1,'Performing Newton method...\n');
end
for i = 1:MAX_ITER
    if i == MAX_ITER, ERROR_FLAG = 1; end
    
    % Form gradient and Hessian
    l = (-A*x+t_min)/sigma; u = (-A*x+t_max)/sigma;
    g = g_quant(l,u)-kappa*g_pen(x);
    H = H_quant(l,u)+kappa*H_pen(x);
    
    % Compute Newton step
    dx = H\(-g);
    
    % Get into feasible region
    gamma = 1;
    if max(x+gamma*dx)>1
        ind_x = find(dx>0);
        gamma = min(gamma,0.99*min((1-x(ind_x))./dx(ind_x)));
    end
    
    if min(x+gamma*dx)<0
        ind_x = find(dx<0);
        gamma = min(gamma,0.99*min(-x(ind_x)./dx(ind_x)));
    end
    
    % Line search
    norm_cur = norm(g);
    while (1)
        l_next = (-A*(x+gamma*dx)+t_min)/sigma; u_next = (-A*(x+gamma*dx)+t_max)/sigma;
        g_next = g_quant(l_next,u_next)-kappa*g_pen(x+gamma*dx);
        norm_next = norm(g_next);
        if norm_next <= (1-alpha*gamma)*norm_cur, break; end
        gamma = beta*gamma;
    end
    
    x=x+gamma*dx;
    
    norm_ratio = norm_next/sqrt(n);
    if verb
        fprintf(1,'Iteration: %2d, Step size: %e, Loss: %e\n',i,gamma,l_y(x))
    end
    
    if norm_ratio<1e-3, break; end
end
if verb, fprintf(1,'Done!\n\n'); end
if ERROR_FLAG==1
    if verb, fprintf(1,'Ran into numerical problems...\n\n'); end
    return
end
end


%% Rounding
[x_sort,ind_x] = sort(x,'descend'); x_cand = []; l_cand = [];
for i = 1:num_sort
    x_cur = zeros(n,1);
    x_cur(ind_x(1:i)) = 1;
    l_cur = l_y(x_cur);
    x_cand = [x_cand x_cur]; l_cand = [l_cand l_cur]; 
end
[l_sort,ind_l] = sort(l_cand,'ascend');
X_amb = x_cand(:,ind_l(1:K)); %get ambiguity set
l_amb = l_sort(1:K);

if verb
    fprintf(1,'After rounding:\n')
    fprintf(1,'Minimum loss in ambiguity set: %e\n',l_amb(1));
    fprintf(1,'Maximum loss in ambiguity set: %e\n\n',l_amb(K));
end


%% Local Optimization
if LO_FLAG
if verb,
    fprintf(1,'Performing local optimization...\n')
end
EXIT_FLAG = 0; iter = 0;
while(~EXIT_FLAG)
    x_cur = X_amb(:,1); x_best = x_cur;
    for i = 1:n
        iter = iter+1;
        x_cur(i) = not(x_cur(i));
        l_cur = l_y(x_cur);
        if any(l_cur<l_amb)
            ind = find(l_cur<l_amb);
            ind = ind(1);
            X_amb = [X_amb(:,1:(ind-1)) x_cur X_amb(:,ind:(end-1))];
            l_amb = [l_amb(:,1:(ind-1)) l_cur l_amb(:,ind:(end-1))];
            if ind==1
                if verb, fprintf(1,'Found new best pattern!\n'); end
            else
                x_cur(i) = not(x_cur(i));
            end
        else
            x_cur(i) = not(x_cur(i));
        end
    end
    if all(x_best == X_amb(:,1)), 
        EXIT_FLAG = 1; 
    end
end
if verb
    fprintf(1,'Perfomed %d iterations. After local optimization:\n',iter);
    fprintf(1,'Minimum loss in ambiguity set: %e\n',l_amb(1));
    fprintf(1,'Maximum loss in ambiguity set: %e\n\n',l_amb(K));
end
end

if verb, fprintf(1,'Finished!\n\n'); end



