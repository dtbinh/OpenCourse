function [X,U] = proj_consistency(X,U,cliques,cliques_index,oindex,ovalues,adddiag)
    % PROJ_CONSISTENCY : Projection onto the consistency set
    %
    %   (X,U) = proj_consistency(X,U) : projects matrix X and vector U onto 
    %   the subspace
    %
    %      {(X,U) : P_k*X*P_k.T = U_k, U = (vec(U_1), ..., vec(U_l)) }.
    %
    % The operation is of the format 
    %
    % INPUT: 
    %   X      	: N x N matrix 
    %   U     	: stacked vector (vec(U_1), ..., vec(U_l)) 
    %   cliques : a cellarray containing the indices beta_k of the cliques
    %   ovalues : a matrix where 
    %             ovalues[i,j] = |{beta_k | i and j in beta_k}|
    %             In other words, the (i,j)th value tells how many times
    %                   element (i,j) appears in a clique
    %   oindex  : index of nonzeros in overlap matrix, in vectorized form
    %             overlap matrix is O = sum_k P_k.T*P_k
    %             then overlap(oindex) = ovalues
    % OUTPUT:
    %   X : matrix output
    %   U : stacked vector (vec(U_1), ..., vec(U_l)) output
    %
    % Authors: Yifan Sun & Lieven Vandenberghe 
    % Date: March 2015

    if size(X,1) > size(X,2)
        N = sqrt(size(X,1));
        X = reshape(X,N,N);
    else
        N = size(X,1);
    end

    if adddiag
        veclen = length(U)-N;
        Xvec = [sparsemat2vec(X,cliques,veclen) ; diag(X)];
        Bvec = Xvec - U;
        Bmat = vec2sparsemat(Bvec,[cliques_index;[[1:N]',[1:N]']],N);
        Bmat(oindex) = Bmat(oindex)./(ovalues+1);
        X = X - Bmat;
        z = Bvec - [sparsemat2vec(Bmat,cliques,veclen); diag(Bmat)];
        U = U + z;
    else
        veclen = length(U);
        Xvec = sparsemat2vec(X,cliques,veclen);
        Bvec = Xvec - U;
        Bmat = vec2sparsemat(Bvec,cliques_index,N);
        Bmat(oindex) = Bmat(oindex)./(ovalues+1);
        X = X - Bmat;
        z = Bvec - sparsemat2vec(Bmat,cliques,veclen);
        U = U + z;
    end

end
