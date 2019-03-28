function val = l2_log(x, auxdata)
    C = auxdata{1};
    z = auxdata{2};
    u = auxdata{3};
    rho = auxdata{4};
    %N = auxdata{5};
    
    %m = size(C,1);
    val = sum(log(1 + exp(C*x))) + (1/2)*(x - z + u)'*rho*(x - z + u);
end