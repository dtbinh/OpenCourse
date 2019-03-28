
function [Sp, Spfill,Sp_um,p, overlap, veclen, cliques, cliques_um] ...
    = chordal_embedding(Spstart,tfill,tsize, display_unmerged)
    % CHORDAL_EMBEDDING : Perform chordal embedding of sparsity pattern
    %
    % 
    %   [Sp, Spfill,Sp_um,p, overlap, veclen, cliques, cliques_um] ...
    %   = chordal_embedding(Spstart,tfill,tsize, display_unmerged):
    %   returns the sparsity pattern of the filled and unfilled patterns,
    %   the cliques (merged and, if desired, unmerged) and other useful
    %   variables.
    %
    % INPUT: 
    %   Spstart          : original sparsity pattern
    %   tfill, tsize     : parameters of merging
    %   display_unmerged : 1 to return unmerged sparsity pattern, 0 else
    %
    % OUTPUT:
    %   Sp      : original sparsity pattern after permutation
    %   Spfill  : sparsity pattern permuted and filled
    %   p       : permutation
    %   overlap : n x n matrix, ij: number of cliques containing i or j
    %   veclen  : sum |beta_k|^2 where beta_k are cliques
    %   cliques : cliques of chordal embedding after permutation
    %
    % Authors: Yifan Sun & Lieven Vandenberghe 
    % Date: March 2015

    [~,~,cliques,~,cliques_um, ~, p,~] = symbolic(Spstart, amd(Spstart), tfill,tsize,display_unmerged);
    
    Sp = Spstart(p,p);
    Spfill = Sp*0;
    for k = 1:length(cliques)
        cl = cliques{k};
        Spfill(cl,cl) = 1;
    end
    if display_unmerged
        Sp_um = Sp*0;
        for k = 1:length(cliques_um)
            cl = cliques_um{k};
            Sp_um(cl,cl) = 1;
        end
    else
        Sp_um = 0;
    end
    
    veclen = 0;
    overlap = Sp*0;
    for k = 1:length(cliques)
        cl = cliques{k};
        ncl = length(cl);
        veclen = veclen + ncl^2;
        overlap(cl,cl) = overlap(cl,cl) + 1;
    end
end
