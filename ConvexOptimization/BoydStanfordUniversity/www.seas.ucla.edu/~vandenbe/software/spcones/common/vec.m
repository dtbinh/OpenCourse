function x = vec(X)
    % VEC : For matrices, returns in vectorized format (stacked columns.)
    %
    %   x = vec(X) takes a matrix X and returns in a vectorized form x
    %
    % INPUT
    %   X : matrix
    %
    % OUTPUT
    %   x : vectorized matrix
    %
    % Authors: Yifan Sun & Lieven Vandenberghe 
    % Date: March 2015
    
    x = X(:);
end
