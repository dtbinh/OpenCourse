function [X, out] = IPM(solver,problemtype, A, D)
    % IPM : Interior point method wrapper (using cvx) for projection onto sparse
    % matrix cones.
    %
    %   [X, out] = IPM((solver,problemtype, A, D): projects sparse matrix 
    %   A onto matrix cones using CVX (interior point solver). Returns X 
    %   the projected matrix and out, a structure containing runtime details.
    %
    % INPUTS
    %   solver      : string indicating which cvx solver to use 
    %                 (e.g. sedumi, sdpt3, mosek). Solver must already be 
    %                 installed.
    %   problemtype : string indicating which problems to solve: 
    %                  'sdp'  : sparse positive semidefinite cone
    %                  'sdpc' : sparse positive semidefintie completable cone
    %                  'edmc' : sparse Euclidean distance matrix
    %                           completable cone
    %   A           : Matrix to project
    %   D           : Sparsity pattern adjacency matrix
    % OUTPUTS
    %   X   : Projection solution
    %   out : Tracked details. Fields include
    %       runtime : CPU time used in solving. Does not include
    %                 preprocessing.
    %       obj     : Final objective value.
    %
    % Authors: Yifan Sun & Lieven Vandenberghe 
    % Date: March 2015
    
    N = size(A,1);
    cvx_solver(solver)
    if strcmpi(problemtype,'sdp')
        cvx_begin 
            variable X(N,N) symmetric
            minimize (1/2*sum_square(vec((X-A).*D)))
            subject to
                X == semidefinite(N)
                X.*D == X
        cvx_end
    elseif strcmpi(problemtype,'sdpc')
        cvx_begin
            variable X(N,N) symmetric
            minimize (1/2*sum_square(vec((X-A).*D)))
            subject to
                X == semidefinite(N)
        cvx_end
    elseif strcmpi(problemtype,'edmc')

        I = [ones(N-1,1) ; [2:N]'];
        J = [[1:N-1]' ; [1:N-1]'];
        V = [-ones(N-1,1) ; ones(N-1,1)];
        V = sparse(I,J,V,N,N-1);

        cvx_begin
            variable X(N,N) symmetric
            minimize (1/2*sum_square(vec((X-A).*D)))
            subject to
                -V'*X*V == semidefinite(N-1)
                diag(X) == 0
        cvx_end
    end
    out.runtime = cvx_cputime;
    out.obj = cvx_optval;
end
