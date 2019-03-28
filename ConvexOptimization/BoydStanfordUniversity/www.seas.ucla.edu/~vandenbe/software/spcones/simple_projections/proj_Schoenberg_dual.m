function X = proj_Schoenberg_dual(X,sp)
    % PROJ_SCHOENBERG_DUAL : Projection onto the dual of dense EDM cone
    %
    %   X = proj_Schoenberg_dual(X, sp) : projects matrix X onto cone 
    %
    %       Z = {X : c'*X*c <= 0 whenever sum(c) = 0} 
    %
    %   This is the dual of the dense matrix cone
    %
    %       S = {VXV' : X <= 0} 
    %
    %   where V is an n x n-1 matrix with orthogonal rows, V*ones(n,1) = 0.
    %   By Schoenberg's condition, X is an EDM iff X is in S and has 0
    %   diagonal.
    %
    %
    % INPUT: 
    %   X   : N x N matrix
    %   sp  : In the eigenvalue decomposition:
    %           1 to use eigs (Lanczos-type)
    %           0 to use eig (full eigenvalue decomposition)
    % OUTPUT:
    %   X   : projected N x N matrix
    %
    % Authors: Yifan Sun & Lieven Vandenberghe 
    % Date: March 2015

    X = .5*(X+X');
    N = size(X,1);
    W = [(-1/sqrt(N)*ones(1,N-1));...
        (speye(N-1)-1/(N+ sqrt(N)))];
    WXW = W'*X*W;
    WXW = .5*(WXW+WXW');
    X = W*proj_nPSD(WXW,sp)*W';
    X = .5*(X+X');
    
end
