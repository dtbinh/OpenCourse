function opts = preprocess_opts(A,opts)
    % PREPROCESS_OPTS : Helper function to format problem data in many 
    % ways, for each solver
    %
    %   opts = preprocess_opts(A,opts) returns options structure with
    %   redundant fields (used for preprocessing) populated
    %
    % INPUT: 
    %   opts: structure containing fields mandatory field:
    %           cliques
    % OUTPUT:
    %   opts: structure containing fields
    %           overlap_index, overlap_values : sparse representation of overlap matrix
    %           overlap_pdiag                 : overlap matrix plus eye (for edm)
    %           overlap_pdiag_values          : sparse values for pdiag
    %           vec_of_cliquelen              : list of lengths of cliques
    %           veclen                        : length of vector representation
    %           diag_index                    : index of diagonal
    %           offset_index                  : list of each subvector
    %                                           index in vector representation
    %           cliques_index                 : sparse index of each submatrix
    %           cliques_vecIJ                 : list of cells containing indices of submatrices
    %           Dindex                        : index of nonzeros in vectorized matrix
    %
    % Authors: Yifan Sun & Lieven Vandenberghe 
    % Date: March 2015


N = size(A,1);
D = opts.D;
A = A.*D;

Dindex = find(sparse(vec(D)));
veclen = 0;
cliques = opts.cliques;
cliquelen = length(cliques);

offset_index = [0];
cliquesI = [];
cliquesJ = [];
overlapIsmall = [];
overlapJsmall = [];
cliques_vecI = [];
cliques_vecJ = [];
cliques_vecIJ = [];
vec_of_cliquelen = zeros(cliquelen,1);
for k = 1:cliquelen
    if mod(k,250) == 0
        cliquesI = [cliquesI ; overlapIsmall];
        cliquesJ = [cliquesJ ; overlapJsmall];
        overlapIsmall = [];
        overlapJsmall = [];
    end
    cl = cliques{k};
    ncl = length(cl);
    vec_of_cliquelen(k) = ncl;
    veclen = veclen + ncl^2;
    indexI =  kron(ones(ncl,1),cl);
    indexJ =  kron(cl,ones(ncl,1));
    overlapIsmall = [overlapIsmall ; indexI];
    overlapJsmall = [overlapJsmall ; indexJ];
    cliques_vecI = [cliques_vecI, {indexI}];
    cliques_vecJ = [cliques_vecJ, {indexJ}];
    cliques_vecIJ = [cliques_vecIJ, {N*(indexJ-1) + indexI}];
    % 		overlap11(cl,cl) = overlap11(cl,cl) + 1;
    offset_index = [offset_index, offset_index(end) + length(cl)^2];
end


cliquesI = [cliquesI ; overlapIsmall];
cliquesJ = [cliquesJ ; overlapJsmall];

overlap = sparse(cliquesI,cliquesJ,1,N,N);

opts.overlap_index = find(sparse(overlap(:)>0));
opts.overlap_values = overlap(opts.overlap_index);

opts.overlap_pdiag = overlap + speye(N);
opts.overlap_pdiag_values = full(opts.overlap_pdiag(opts.overlap_index));

opts.vec_of_cliquelen = vec_of_cliquelen;
opts.veclen = veclen;
opts.diag_index = find(speye(N));
opts.offset_index = offset_index;
opts.cliques_index = [cliquesI,cliquesJ];
opts.cliques_vecIJ = cliques_vecIJ;
opts.Dindex = Dindex;
end
