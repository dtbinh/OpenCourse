function y = sparsemat2vec(X,cliques,veclen)

    % SPARSEMAT2VEC : For matrix cones, selects lower dimensional cones.
    %
    %   y = sparsemat2vec(X,cliques,veclen) : Takes X an N x N matrix, and 
    %   returns y = (y_1,...,y_l) the vectorized versions of y_k, where,
    %   for cliques bk, k = 1,..., l,
    %       y_k = vec(P_k*X*P_k.T),  P_k X P_k.T = X(bk,bk).
    %
    %INPUT: 
    %   X       : N x N sparse matrix input
    %   cliques : a cellarray containing the indices beta_k of the cliques
    %   veclen  : length of y, = sum_k |bk|
    %
    % OUTPUT:
    %   y       : vectorized output (y_1,...,y_l) where y_i are the vector 
    %             formats of the yth lower dimensional cone variable
    %
    % Authors: Yifan Sun & Lieven Vandenberghe 
    % Date: March 2015

    y = zeros(veclen,1);
    if size(X,1) > size(X,2)
        N = sqrt(size(X,1));
        X = reshape(X,N,N);
    end
    offset = 0;
    for k = 1:length(cliques)
        cl = cliques{k};
        ncl = length(cl);
        y(offset +1: offset + ncl^2) = vec(X(cl,cl));
        offset = offset + ncl^2;
    end
end
