function [time,obj,itr] = solve_single_prox_via_cvx(A,b,C,Z, solver)
    
    % Solves
    % 
    %     minimize trace(C*X) + 1/2 * ||X-Z||_F^2
    %     subject to
    %         A(X) == b
    %         X >= 0
    %
    % using CVX
    
    p = size(C,1);
    E = ones(p);
    I = find(tril(E, -1));  


    cvx_solver(solver)
        cvx_begin 
        variable X(p,p) symmetric
        
        minimize (trace(C*X) + 0.5*sum_square([diag(X) ; X(I)*sqrt(2)]-[diag(Z) ; Z(I)*sqrt(2)]))
        subject to
            A'*X(:) == b
            X == semidefinite(p)

    cvx_end
    
    time = cvx_cputime;
    obj = cvx_optval;
    itr = cvx_slvitr;
  
end
