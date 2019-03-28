function [Xmat,out] = proxgrad(A,opts)
    % PROXGRAD : Proximal gradient solver to project onto the sparse 
    % positive semidefinite, positive semidefinite completable, or 
    % Euclidean distance completable matrices.
    %
    %   [X, out] = proxgrad(A, opts) : projects sparse matrix A onto matrix
    %   cones using the proximal gradient method, with details embedded in 
    %   structure opts. Returns X the projected matrix and out, a structure 
    %   containing runtime details.
    %
    %
    % INPUTS
    %   A       : N x N sparse matrix to be projected
    %   opts    : structure with input options. Required fields inputed by
    %             user include
    %       problemtype : 'sdp' if projecting onto the sparse positive
    %                     semidefinite cone, 
    %                     'sdpc' if projecting onto the cone of matrices
    %                     with a positive semidefinite completion
    %                     'edm' if projecting onto the cone of matrices
    %                     with an EDM completion
    %       D           : adjacency (binary) of sparsity pattern
    %       cliques     : a cellarray containing the indices beta_k of the
    %                     cliques of E (if E is chordal) or the embedding of E
    %                     (if E is nonchordal).
    %       accelerate  : boolean, 1 to use fast prox gradient, 
    %                              0 to use default prox grad
    %       tol         : tolerance for stopping condition
    %       maxiter     : maximum number of iterations
    %       maxtime     : maximum seconds allowed
    %       verbose     : boolean, 1 for printout while running, 0 else
    %       epoch       : number of iterations between checking stopping condition
    %
    %             Additional fields supplied by preprocessing (user can
    %             ignore these, just use preprocess_opts)
    %       offset_index     : list containing pointers for where each
    %                          subvector Zk appears in Z, where Zk corresponds to submatrix
    %                          X_{bk,bk}
    %       vec_of_cliquelen : vector containing size of each clique
    %       cliques_index    : m x 2 matrix pertaining to the row and
    %                          column indices corresponding to each
    %                          submatrix X_{bk,bk}. 
    %       veclen           : length of vector Z
    %       L                : Lipschitz constant of objective, equal to
    %                          max number of overlaps of indices in cliques
    %       
    %       
    % OUTPUTS
    %   X    : N x N solution to matrix nearness problem
    %   out  : structure with tracking details, returned by solver
    %          obj     : primal objective, calculated at every iteration
    %          err     : residuals, calculated at every epoch.
    %          maxiter : total number of iterations before termination.
    %          runtime : total cpu time needed for internal solver
    %
    % Authors: Yifan Sun & Lieven Vandenberghe 
    % Date: March 2015
    
    % Details on algorithm: The problem solved is 
    %
    %       minimize ||sum_k Pk.T*Xk*Pk-A||_F^2
    %       subject to Xk \in C_k
    %
    % where Ck is either a PSD or EDM cone of order |bk|.
    %
    % We use this method to solve three problems:
    %
    %   1) projection onto the sparse SDP cone (sdp)
    %
    %   2) projection onto the cone of sparse matrices with SDP completion (sdpc)
    %
    %   3) projection onto the cone of sparse matrices with EDM completion (edm)
    %
    % where the sparsity pattern E is chordal with cliques bk. 
    %
    % For problem 1: the method solves the stated problem, with Ck = PSD cones.
    %
    % For problems 2 and 3: the method solves the dual problem (dualmode),
    % with Ck = PSD or EDM cones (resp)
    


accelerate = opts.accelerate;
tol = opts.tol;
maxiter = opts.maxiter;
verbose = opts.verbose;
epoch = opts.epoch;

offset_index = opts.offset_index;
vec_of_cliquelen = opts.vec_of_cliquelen;
cliques_index = opts.cliques_index;
cliques = opts.cliques;
veclen = opts.veclen;
cliqueslen = length(cliques);
N = size(A,1);


L = opts.L;
if strcmpi(opts.problemtype,'sdp') ||  strcmpi(opts.problemtype,'sdpc')
    X = sparsemat2vec(A,cliques,N);
    ptype = 1;
elseif strcmpi(opts.problemtype,'edmc')
    cliques_index = [cliques_index;[[1:N]',[1:N]']];
    X = [sparsemat2vec(A,cliques,N) ; zeros(N,1)];
    L = L + 1;
    ptype = 2;
end

if strcmpi(opts.problemtype,'edmc') ||  strcmpi(opts.problemtype,'sdpc')
    dualmode = true;
    A = -A;
elseif strcmpi(opts.problemtype,'sdp')
    dualmode = false;
end

t = 1/L;

obj = zeros(maxiter,1);
D = opts.D;
if accelerate
    V = X;
    Xp = X;
else
    Xp = X;
end

if verbose
    fprintf('\n\n')
    if accelerate
        fprintf('Fast proximal ');
    else
        fprintf('Proximal ');
    end
    fprintf('gradient method to solve the projection onto\n');
    
    if strcmpi(opts.problemtype,'edmc')
        fprintf('the cone of sparse matrices with EDM completion.\n');
    elseif strcmpi(opts.problemtype,'sdpc')
        fprintf('the cone of sparse matrices with PSD completion.\n');
    elseif strcmpi(opts.problemtype,'sdp')
        fprintf('the cone of sparse PSD matrices.\n');
    end
    if dualmode
        fprintf('The method is solving the dual problem.\n');
    end
    fprintf('\n')

end
starttime = cputime;
for iter = 1:maxiter
    
    %%smoothF
    if accelerate
        theta = 2/(iter+1);
        Y = (1-theta)*X+theta*V;
        smoothterm = Y;
    else
        smoothterm = X;
    end
    
    xmat = vec2sparsemat(smoothterm,cliques_index,N);
    
    if ptype == 1
        dgmat = xmat-A;
        dg = sparsemat2vec(dgmat,cliques,veclen);
    elseif ptype == 2
        dgmat = xmat-A;
        dg = full([sparsemat2vec(dgmat,cliques,veclen);diag(dgmat)]);
    end
    
    proxterm = smoothterm-dg*t;
    
    %%project step
    for k = 1:cliqueslen
        ncl = vec_of_cliquelen(k);
        index = (offset_index(k) + 1):(offset_index(k) + ncl^2);
        if ptype == 1
            tmp = proj_PSD(reshape(proxterm(index),ncl,ncl),0);
        elseif ptype == 2
            tmp = proj_Schoenberg_dual(reshape(proxterm(index),ncl,ncl),0);
        end
        X(index) = tmp;
    end
    if ptype == 2
        X(veclen+1:end) = proxterm(veclen+1:end);
    end
    
    
  
    % check stopping condition     
    if mod(iter,epoch) == 0 || (epoch == 1)
        Xmat = vec2sparsemat(X,cliques_index,N);
        if dualmode
            obj(iter) = sum(sum((D.*Xmat).^2))/2;
        else
            obj(iter) = sum(sum((D.*(Xmat-A)).^2))/2;
        end
        
        
        if accelerate
            r = (Y-X)/t;
        else
            r = (X-Xp)/t;
        end
        err(iter).primal = norm(r,'fro')/max(norm(X,'fro'),1);
        
        
        if verbose
            fprintf('it:%d, obj: %.4e, ep: %.4e, \n',iter,obj(iter), err(iter).primal);            
        end
        endtime = (cputime - starttime);
        if ((err(iter).primal<tol) && (iter >= epoch)) || (endtime > opts.maxtime)
            err = err(1:iter);
            obj = obj(1:iter);
            out.runtime = endtime;
            break
        end
    end
    
    if accelerate
        V = X + 1/theta*(X-Xp);
    end
    Xp = X;
    
    
end
out.obj = obj;
out.err = err;
out.maxiter = iter;

if dualmode
    Xmat = Xmat - A;
end
if strcmpi(opts.problemtype,'edmc')
    Xmat = Xmat - spdiag(diag(Xmat));
end
end
