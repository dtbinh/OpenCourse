function X = proj_nPSD(X,sp)
    % PROJ_NPSD : Projection onto the dense negative semidefinite cone
    %
    %   X = proj_nPSD(X, sp) : projects matrix X onto dense negative
    %   semidefinite cone.
    %
    % INPUT: 
    %   X   : N x N matrix
    %   sp  : In the eigenvalue decomposition, 1 to use eigs (Lanczos-type), 
    %         0 to use eig (full eigenvalue decomposition)
    % OUTPUT:
    %   X   : projected N x N matrix
    %
    % Authors: Yifan Sun & Lieven Vandenberghe 
    % Date: March 2015

    X = proj_PSD(-X,sp);
    X = -X;
end
