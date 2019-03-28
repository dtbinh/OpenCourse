function y = spdiag(x)
    %SPDIAG : Input vector, return sparse matrix with vector along diagonal
    %
    % y = spdiag(x) returns a sparse matrix y with x along its diagonal
    %
    % INPUT
    %   x : vector of length n
    %
    % OUTPUT
    %   y : sparse n x n matrix where diag(y) = x
    %
    % Authors: Yifan Sun & Lieven Vandenberghe 
    % Date: March 2015
    
    n = length(x);
    y= sparse(1:n,1:n,x,n,n);
end
