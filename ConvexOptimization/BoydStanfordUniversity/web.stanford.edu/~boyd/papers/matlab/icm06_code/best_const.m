function [w] = best_const(A)
% BEST_CONST	Returns the constant edge weight that yields fastest 
%               averaging on a graph.
%
% For more details see the references:
% "Fast linear iterations for distributed averaging" by L. Xiao and S. Boyd
% "Convex Optimization of Graph Laplacian Eigenvalues" by S. Boyd
%
% [w] = BEST_CONST(A) gives a vector of the best constant edge weights
% for a graph described by the incidence matrix A (n x m).
% Here n is the number of nodes and m is the number of edges in the graph;
% each column of A has exactly one +1 and one -1.
%
% The best constant edge weight is the inverse of the average of
% the second smallest and largest eigenvalues of the unweighted Laplacian:
%
%   w = 2/( lambda_2(A*A') + lambda_n(A*A') )
%
% Almir Mutapcic 08/29/06
[n,m] = size(A);

% max degrees of the nodes
Lunw = A*A';                % unweighted Laplacian matrix
eigvals = sort(eig(Lunw));

% max degree weigth allocation
alpha = 2/(eigvals(2) + eigvals(n));
w = alpha*ones(m,1);
