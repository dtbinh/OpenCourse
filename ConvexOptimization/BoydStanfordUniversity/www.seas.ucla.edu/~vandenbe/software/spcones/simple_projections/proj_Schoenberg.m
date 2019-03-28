function X = proj_Schoenberg(X,sp)
    % PROJ_SCHOENBERG : Projection onto the dense EDM cone
    %
    %   X = proj_Schoenberg(X, sp) : projects matrix X onto dense matrix cone
    %
    %       S = {X : c'*X*c <= 0 whenever sum(c) = 0}
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

  
    N = size(X,1);
    if N == 1
        X = 0;
        return;
    end
    
    X = .5*(X+X');
    W = [(-1/sqrt(N)*ones(1,N-1));...
        (eye(N-1) - 1/(N+ sqrt(N)))];
    X = X - W*proj_PSD(W'*X*W,sp)*W';
    X = .5*(X+X');
    
end
