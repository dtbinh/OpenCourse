function [w] = fmmc(A)
% FMMC    Returns the fastest mixing Markov chain (FMMC) edge weights
%         (i.e., edge transition probabilities).
%
% [w] = FMMC(A) gives a vector of the fastest mixing Markov chain
% edge weights for a graph described by the incidence matrix A (n x m).
% Here n is the number of nodes and m is the number of edges in the graph;
% each column of A has exactly one +1 and one -1.
%
% The FMMC edge weights are given the SDP:
%
%   minimize    s
%   subject to  -s*I <= I - L - (1/n)11' <= s*I
%               w >= 0,  diag(L) <= 1
%
% where the variables are edge weights w in R^m and s in R.
% Here L is the weighted Laplacian defined by L = A*diag(w)*A'.
%
% For more details see references:
% "Fastest mixing Markov chain on a graph" by S. Boyd, P. Diaconis, and L. Xiao
% "Convex Optimization of Graph Laplacian Eigenvalues" by S. Boyd
%
% Written for CVX by Almir Mutapcic 08/29/06
[n,m] = size(A);

cvx_begin
  variable w(m,1)   % edge weights
  variable s        % epigraph variable

  minimize( s )
  subject to
    -eye(n) + A*diag(w)*A' + ones(n,n)/n + s*eye(n) == semidefinite(n);
    eye(n) -  A*diag(w)*A' - ones(n,n)/n + s*eye(n) == semidefinite(n);

    w >= 0;
    diag(A*diag(w)*A') <= 1;
cvx_end
