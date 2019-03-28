function [X,out] = dr(A,opts)
    % DR : Douglas Rachford solver for finding the projection onto the sparse 
    % positive semidefinite, positive semidefinite completable, or 
    % Euclidean distance completable matrices.
    %
    %   [X, out] = dr(A, opts) : projects sparse matrix A onto matrix
    %   cones using the Douglas-Rachford method, with details embedded in 
    %   structure opts. Returns X the projected matrix and out, a structure 
    %   containing runtime details.
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
    %       Dindex           : index of nonzeros in vectorized sparse matrix
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
            
    
    % Details on algorithm: The primal problem solved is 
    %
    %   A)    minimize sum_{ij\in E} (X_{ij}-A_{ij})^2
    %         subject to X_{bk,bk} \in C_k
    %
    % where Ck is either a PSD or EDM cone of order |bk| and bk are the 
    % cliques of the sparisty pattern. 
    %
    % The dual problem can be written as
    %
    %   B)    minimize ||X+A||_F^2
    %         subject to Zk \in Ck
    %                    X = sum_k Pk.T * Zk* Pk
    %                    X has sparsity E
    %
    % We use this method to solve three problems:
    %
    %   1) projection onto the sparse SDP cone (sdp)
    %
    %       minimize    ||X-A||_F^2
    %       subject to  X >= 0,
    %                   X in sparsity pattern E'
    %
    %   2) projection onto the cone of sparse matrices with SDP completion (sdpc)
    %
    %       minimize ||X-A||_F^2
    %       subject to X in sparsity pattern E,'
    %                  X has a PSD completion
    %
    %   3) projection onto the cone of sparse matrices with EDM completion (edm)
    %
    %       minimize ||X-A||_F^2
    %       subject to X in sparsity pattern E',
    %                  X has an EDM completion
    %
    % where the sparsity pattern E has chordal extension E' with cliques bk. 
    %
    % For problem 1, we solve problem B or problem A in dual mode. For
    % problems 2 and 3, we solve problem A or problem B in dual mode.


t = opts.t;
rho = opts.rho;

maxiter = opts.maxiter;
verbose = opts.verbose;
epoch = opts.epoch;

tol = opts.tol;
veclen = opts.veclen;
cliques = opts.cliques;
cliqueslen = length(cliques);
cliques_index = opts.cliques_index;
overlap_index = opts.overlap_index;
overlap_values  = opts.overlap_values;
D = opts.D;
offset_index = opts.offset_index;
vec_of_cliquelen = opts.vec_of_cliquelen;
dualmode = opts.isdual;
Dindex = opts.Dindex;
obj = zeros(maxiter,1);


Xvec = zeros(veclen,1);
Zvec = zeros(veclen,1);

N = size(A,1);
if strcmpi(opts.problemtype,'sdp') 
    ptype = 1;
elseif strcmpi(opts.problemtype,'sdpc')
    ptype = 2;
elseif strcmpi(opts.problemtype,'edmc')
    ptype = 3;
    overlap_pdiag_values = opts.overlap_pdiag_values;
    diag_index = opts.diag_index;
    Xvec = [Xvec ; zeros(N,1)];
    Zvec = [Zvec ; zeros(N,1)];
end

if dualmode
    A = -A;
end


Xmat = A;
Zmat = A;


if verbose

    fprintf('\n\n')
    fprintf('Douglas Rachford method to solve the projection onto\n');
    
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
    %X = proxA(Z,t);
    %update xvec
    for k = 1:cliqueslen
        ncl = vec_of_cliquelen(k);
        toproject = reshape(Zvec(offset_index(k) + 1:offset_index(k) + ncl^2),ncl,ncl);
        if ptype==1 || ptype == 2
            tmp = proj_PSD(toproject,0);
        elseif ptype == 3 && ~dualmode
            tmp = proj_Schoenberg(toproject,0);
        elseif ptype == 3 && dualmode
            tmp = proj_Schoenberg_dual(toproject,0);
        end
        tmp = vec(tmp);
        Xvec(offset_index(k) + 1:offset_index(k) + ncl^2) = tmp;
    end

       
    %update xmat
    rhs = (D.*A+t*Zmat)/(1+t);
    if (ptype == 1 && ~dualmode) || (ptype == 2 && dualmode)
        Xmat = rhs;
    elseif (ptype == 2 && ~dualmode) || (ptype == 1 && dualmode)
        Xmat = Zmat;
        Xmat(Dindex) = rhs(Dindex);
    elseif  ptype == 3 && ~dualmode
        Xvec(veclen+1:end) = 0;
        Xmat = Zmat;
        Xmat(Dindex) = rhs(Dindex);
    elseif ptype == 3 && dualmode
        Xmat = rhs;
        Xvec(veclen+1:end) = Zvec(veclen+1:end);
    end
    

    if mod(iter,epoch) == 0
        XZmat = Xmat-Zmat;
        XZvec = Xvec-Zvec;
        
        if (ptype == 1 && ~dualmode) || (ptype == 2 && dualmode)
            [pXmat,pXvec] = proj_sum_consistency(Xmat,Xvec,cliques,cliques_index,overlap_index,overlap_values,Dindex);
        elseif (ptype == 2 && ~dualmode) || (ptype == 1 && dualmode)
            [pXmat,pXvec] =  proj_consistency(Xmat,Xvec,cliques,cliques_index,overlap_index,overlap_values,0);
        elseif ptype == 3 && dualmode
            [pXmat,pXvec] =  proj_sum_consistency(Xmat,Xvec,cliques,cliques_index,overlap_index,overlap_values,Dindex, diag_index);
        elseif ptype == 3 && ~dualmode
            [pXmat,pXvec] =  proj_consistency(Xmat,Xvec,cliques,cliques_index,overlap_index,overlap_pdiag_values,1);
        end
    end
    
    
    Ymat = 2*Xmat - Zmat;
    Yvec = 2*Xvec - Zvec;
    if (ptype == 1 && ~dualmode) || (ptype == 2 && dualmode)
        [Ymat,Yvec] = proj_sum_consistency(Ymat,Yvec,cliques,cliques_index,overlap_index,overlap_values,Dindex);
    elseif (ptype == 2 && ~dualmode) || (ptype == 1 && dualmode)
        [Ymat,Yvec] =  proj_consistency(Ymat,Yvec,cliques,cliques_index,overlap_index,overlap_values,0);
    elseif ptype == 3 && dualmode
        [Ymat,Yvec] =  proj_sum_consistency(Ymat,Yvec,cliques,cliques_index,overlap_index,overlap_values,Dindex, diag_index);
    elseif ptype == 3 && ~dualmode
        [Ymat,Yvec] =  proj_consistency(Ymat,Yvec,cliques,cliques_index,overlap_index,overlap_pdiag_values,1);
    end
    
    
    
    Zmat = Zmat + rho*(Ymat-Xmat);
    Zvec = Zvec + rho*(Yvec-Xvec);
    
    
    % update
    if mod(iter,epoch) == 0
        
        
        if dualmode
            obj(iter) = sum(sum((D.*(Xmat)).^2))/2;
        else
            obj(iter) = sum(sum((D.*(pXmat-A)).^2))/2 ;
        end
        
        
        
        
        nX = sqrt(max(norm(Xmat,'fro')^2 + norm(Xvec,'fro')^2,1));
        nXZ = sqrt(max(norm(XZmat,'fro')^2 + norm(XZvec,'fro')^2,1));
        errprim = sqrt(norm(Xmat-pXmat,'fro')^2 + norm(Xvec-pXvec,'fro')^2);
        
        errdual = sqrt(norm(Ymat-pXmat,'fro')^2 + norm(Yvec-pXvec,'fro')^2);
        err(iter).primal = errprim/nX;
        err(iter).dual = errdual/nXZ;
        
        
        if verbose
            fprintf('it:%d, obj: %.4e, ep: %.4e, ed: %.4e \n',iter,obj(iter),full(err(iter).primal),full(err(iter).dual));
        end
        endtime = (cputime - starttime);
        out.runtime(iter) = endtime;
        if((err(iter).primal < tol)  &&(err(iter).dual < tol)&& iter >= epoch) || (endtime > opts.maxtime)
            obj = obj(1:iter);
            err = err(1:iter);
            break
        end
    end
    
    
end
out.obj = obj;
out.err = err;
X = pXmat;
if dualmode
    X = D.*(X - A);
end
out.maxiter = iter;
end
