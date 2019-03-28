function [X,U] = proj_sum_consistency(X, U, cliques, cliques_index, oindex, ovalues, sigma,eta)
    % PROJ_SUM_CONSISTENCY : Projection onto the dual of the consistency set.
    %
    %   (X,U) = proj_sum_consistency(X,U) : projects matrix X and vector U 
    %   onto the subspace
    %
    %     V = {(X,U) :{(X,u) : sum_k P_k.T*mat(u_k)*P_k = X, U = (u_1, ..., u_l) }
    %
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

if exist('sigma','var')
    S = reshape(sparse(sigma,1,1,N^2,1),N,N);
    X = sparse(X.*S);
end

if exist('eta','var')
    veclen = length(U) - length(eta);
    V = U(veclen+1:end);
    U = U(1:veclen);
else
    veclen = length(U);
end

Wmat = sparse(X - vec2sparsemat(U,cliques_index,N));
            
            
if exist('eta','var')
    Wmat(eta) = Wmat(eta) - V;
end

if exist('sigma','var')
	oindex2 = [oindex ; sigma];
	ovalues2 = [ovalues ; ones(size(sigma))];
	
    %overlap(sigma) = overlap(sigma) + 1;
	
    if exist('eta','var')
		oindex2 = [oindex2 ; eta];
		ovalues2 = [ovalues2 ; ones(size(eta))];
        %overlap(eta) = overlap(eta) + 1;
    end
	overlap = reshape(sparse(oindex2,1,ovalues2,N^2,1),N,N);
	Wmat(oindex) = Wmat(oindex)./(overlap(oindex));
	X = X - Wmat.*S;
	
else
    Wmat(oindex) = Wmat(oindex)./(1+ovalues);
    X = X - Wmat;
end

U = U + sparsemat2vec(Wmat,cliques,veclen);

if exist('eta','var')
    V = V + Wmat(eta);
    U = [U ; V];
end

end
