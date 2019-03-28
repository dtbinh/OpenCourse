function Y = vec2sparsemat(x,cliques_index,N)

    % VEC2SPARSEMAT : For matrix cones, expands and adds lower dimensional 
    % cones.
    %
    %   Y = vec2sparsemat(x,cliques_index,N) : Takes x = (x_1,...,x_l)  
    %   the vectorized versions of the lower dimensional conic variables 
    %   and returns Y the N x N sparse matrix, where.
    %
    %     Y = sum_k P_k.T*mat(x_k)*P_k 
    %
    % where P_k*Y*P_k.T selects submatrix Y(beta_k,beta_k), and 
    % P_k.T*Y*P_k is the adjoint operator
    %
    %INPUT: 
    %   x       : vectorized input of (x_1,...,x_l) where x_i are the vector 
    %             formats of the ith lower dimensional cone variable
    %   cliques : a cellarray containing the indices beta_k of the cliques
    %   N       : the order of the matrix Y
    %
    % OUTPUT:
    %   Y       : the solution as an N x N sparse matrix
    %
    % Authors: Yifan Sun & Lieven Vandenberghe 
    % Date: March 2015
    
   Y =  sparse(cliques_index(:,1), cliques_index(:,2),x,N,N);
end
