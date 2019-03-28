function [X,out] = bcd(A,opts)

    % BCD : Block coordinate descent solver for finding the projection onto the 
    % sparse positive semidefinite, positive semidefinite completable, or 
    % Euclidean distance completable matrices.
    %
    %   [X, out] = bcd(A, opts) : projects sparse matrix A onto matrix
    %   cones using the block coordinate descent method (Dykstra's method), 
    %   with details embedded in structure opts. Returns X the projected 
    %   matrix and out, a structure containing runtime details.
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
    %       cliques_vecIJ    : list of cells of I and J indices for each
    %                          clique
    %       veclen           : length of vector Z
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
    
    % The problem solved is 
    %
    %       minimize ||X-A||_F^2
    %       subject to X_{bk,bk} \in C_k
    %
    % where Ck is either a PSD or EDM cone of order |bk| and bk are the 
    % cliques of the sparisty pattern.
    %
    % We use this method to solve three problems:
    %
    %   1) projection onto the sparse SDP cone (sdp)
    %
    %       minimize    ||X-A||_F^2
    %       subject to  X >= 0,
    %                   X in sparsity pattern E
    %
    %   2) projection onto the cone of sparse matrices with SDP completion (sdpc)
    %
    %       minimize ||X-A||_F^2
    %       subject to X in sparsity pattern E,
    %                  X has a PSD completion
    %
    %   3) projection onto the cone of sparse matrices with EDM completion (edm)
    %
    %       minimize ||X-A||_F^2
    %       subject to X in sparsity pattern E,
    %                  X has an EDM completion
    %
    % where the sparsity pattern E is chordal with cliques bk. 
    %
    % For problem 1: the method solves the dual problem (dualmode), 
    % with Ck = PSD cones.
    %
    % For problems 2 and 3: the method solves the primal problem 
    % with Ck = PSD or EDM cones (resp)

maxiter = opts.maxiter;
epoch = opts.epoch;
verbose = opts.verbose;
obj = zeros(maxiter,1);
offset_index = opts.offset_index;
D = opts.D;
cliques = opts.cliques;
cliquelen = length(cliques);
cliques_vecIJ = opts.cliques_vecIJ;
veclen = opts.veclen;
N = size(A,1);
if strcmpi(opts.problemtype,'sdp') || strcmpi(opts.problemtype,'sdpc')
    ptype = 1;
    ncones = cliquelen;
    Z = zeros(veclen,1);
elseif strcmpi(opts.problemtype,'edmc')
diag_index = opts.diag_index;
    ptype = 2;
    ncones = cliquelen + 1;
    Z = zeros(veclen+N,1);
end

if  strcmpi(opts.problemtype,'sdpc') || strcmpi(opts.problemtype,'edmc')
    dualmode = false;
elseif strcmpi(opts.problemtype,'sdp')
    dualmode = true;
    A = -A;
end

X = sparse(A);

if verbose

    fprintf('\n\n')
    fprintf('Block coordinate descent method to solve the projection onto\n');
    
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
start = cputime;
for iter = 1:maxiter
    
    Xdiff = 0;
    zfill = [];
    
    startindex = 1;
    for k = 1:ncones
        %update X
        if k <= cliquelen
            cl = cliques{k};
            ncl = length(cl);
           
            clvec = cliques_vecIJ{k};
            xcurr = full(X(clvec));
           
            % use nan instead of 0 when projecting diagonal since otherwise
            % sparse indices need to be reallocated at each iteration.
            if ptype == 2
                xcurr(isnan(xcurr)) = 0;
            end
          
          
            xv = reshape(xcurr - Z(offset_index(k) + 1 : offset_index(k) + ncl^2),ncl,ncl);
            if ptype == 1
                xv = proj_PSD(xv,0);
            elseif ptype == 2
                xv = proj_Schoenberg(xv,0);
            end
            
            Xd = vec(xv) - xcurr;
            if ncl == 1 && ptype == 2
                X(clvec) = nan;
            else
                X(clvec) = xv;
            end
            
        else
            Xd = -X(diag_index);
            X = X + spdiag(nan*ones(N,1));
        end
        
        
        % update Z
        if  (mod(iter,epoch) == 0 )
            Xdiff = Xdiff + sum(Xd(:).^2);
        end
        
        if length(zfill) > 1000
            Z(offset_index(startindex) + 1 : offset_index(k)) = Z(offset_index(startindex) + 1 : offset_index(k))  + zfill;
            zfill = [];
            startindex = k;
        end
        zfill = [zfill ; Xd(:)];
    end
    Z(offset_index(startindex) + 1 : end) = Z(offset_index(startindex) + 1 : end)  + zfill;
    
 
    endtime = cputime - start;
    out.runtime(iter) = endtime;
    if  (mod(iter,epoch) == 0 ) 
        if dualmode
            obj(iter) = sum(sum((D.*X).^2))/2;
        else
            if ptype == 2
                XA = X; 
                XA(isnan(XA)) = 0; 
                obj(iter) = sum(sum((D.*(XA-A)).^2))/2;
                clear XA
            else
                obj(iter) = sum(sum((D.*(X-A)).^2))/2;
            end
        end
        err(iter).primal = sqrt(Xdiff)/max(norm(Z,'fro'),1);
        if verbose
            fprintf('it: %d: obj %f, err %d\n',iter,obj(iter), err(iter).primal);
        end
        if (maxiter > (epoch-1) && err(iter).primal < opts.tol) || (endtime > opts.maxtime)
            obj = obj(1:iter);
            err = err(1:iter);
            break
        end
    end
end

out.obj = obj;
out.err = err;
out.maxiter = iter;
X(isnan(X)) = 0;

if dualmode
    X = X - A;
end
end
