function grad = l2_log_grad(x, auxdata)
    C = auxdata{1};
    z = auxdata{2};
    u = auxdata{3};
    rho = auxdata{4};
    %N = auxdata{5};
    
    %m = size(C,1);
    e2Cx = exp(C*x);
    grad = (C'*(e2Cx./(1 + e2Cx))) + rho*(x - z + u); %1/(m*N)* (C'*(e2Cx./(1 + e2Cx))) + rho*(x - z + u);
end