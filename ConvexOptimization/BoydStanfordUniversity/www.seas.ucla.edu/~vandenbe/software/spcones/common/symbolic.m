function [snpost,snpar,cliques,supernodes,cliques_nm, supernodes_nm, p,ip] = symbolic(A, p, tfill,tsize,display_unmerged)
    
%   SYMBOLIC: performs symbolic factorization on sparse matrix.
%
%       [snpost,snpar,cliques,supernodes,cliques_nm, supernodes_nm, p,ip] = 
%       symbolic(A, p, tfill,tsize,display_unmerged) : Takes sparse matrix,
%       fill-reducing permutation, and other input parameters, and returns
%       symbolic factorization variables.
%
%   INPUTS:
%       A                : sparse matrix
%       p :                fill-reducing permutation
%       tfill,tsize      : parameters for clique merging. Set to 0 for no
%                          merging.
%       display_unmerged : boolean (1 to return unmerged clique, 0
%                          otherwise.) Depending on clique statistics, 
%                          setting to 1 may significantly increase 
%                          computational complexity.
%
%   OUTPUTS:
%       snpost        : postordering of cliques
%       snpar         : parent of cliques (0 for root)
%       cliques       : cliques (merged)
%       supernodes    : supernodes = setminus(clique, parent clique)
%       cliques_nm    : unmerged cliques
%       supernodes_nm : unmerged supernodes
%       p             : fill reducing permutation
%       ip            : inverse permutation
%
% This script and all those in the "symbolic" directory was adapted from 
% python code in the CHOMPACK package (http://chompack.readthedocs.org/en/latest/)
%

    
    

    [m,n] = size(A);
    if m ~= n
        error('matrix A must be square')
    end

    if isempty(p)
        p = 1:n;
    end
    if nargin == 2
        tfill = 0; tsize = 0;
    end
    
   %% init

    % Symmetrize A
    Ap = tril(A) + tril(A,-1)';

    % Permute if permutation vector p or ordering routine is specified
    if length(p) ~= n
        error('length of permutation vector must be equal to the order of A')
    end
    Ap = A(p,p);
    
    ip = 0*p;
    for k = 1:n
        ip(p(k)) = k;
    end
    
    
    
    [rowidx,colptr,val,n,m] = sparse2ccs(Ap);
    clear Ap
    Ap.rowidx = rowidx;
    Ap.colptr = colptr;
    Ap.val = val;
    Ap.n = n; Ap.m = m;
    

    % Symbolic factorization
    par = myetree(Ap);
    post = post_order(par);
    colcount = counts(Ap, par, post);
    

    [snode, snptr, snpar] = get_supernodes(par, post, colcount);
    snpost = post_order(snpar);

    if display_unmerged
        [sncolptr_nm, snrowidx_nm] = embed(Ap, colcount, snode, snptr, snpar, snpost);
        for k = 1:length(snpar)
            cliques_nm(k) = {snrowidx_nm(sncolptr_nm(k):sncolptr_nm(k+1)-1)};
            supernodes_nm(k) = {snode(snptr(k):snptr(k+1)-1)};
        end
    else
        cliques_nm = []; supernodes_nm = [];
    end
    if (tfill > 0) || (tsize > 0)
        [colcount, snode, snptr, snpar, snpost] = amalgamate(colcount, snode, snptr, snpar, snpost, tfill, tsize);
    end

    [sncolptr, snrowidx] = embed(Ap, colcount, snode, snptr, snpar, snpost);
%     [relptr, relidx] = relative_idx(sncolptr, snrowidx, snptr, snpar);
%     for k = 1:length(snpar)
%      
%         eval(sprintf('cliques%d = snrowidx(sncolptr(k):sncolptr(k+1)-1);',k-1))
%         eval(sprintf('supernodes%d = snode(snptr(k):snptr(k+1)-1);',k-1))
%     end
% save('current2.mat','sncolptr','snrowidx','snpost','snode','snptr','snpar','par','post','colcount','rowidx','colptr','val','p','ip','m','n','cliques*','supernodes*','p','ip')
    for k = 1:length(snpar)
        cliques(k) = {snrowidx(sncolptr(k):sncolptr(k+1)-1)};
        supernodes(k) = {snode(snptr(k):snptr(k+1)-1)};

    end
end
   
    
