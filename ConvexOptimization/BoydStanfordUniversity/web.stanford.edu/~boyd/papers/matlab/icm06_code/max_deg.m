function [w] = max_deg(A)
% MAX_DEG    Returns the maximum-degree heuristic edge weights.
%
% [w] = MAX_DEG(A) gives a vector of maximum-degree edge weights
% for a graph described by the incidence matrix A (n x m).
% Here n is the number of nodes and m is the number of edges in the graph;
% each column of A has exactly one +1 and one -1.
%
% Maximum-degree edge weights are all equal to one over the maximum
% degree of the nodes in the graph.
%
% For more details see the references:
% "Fast linear iterations for distributed averaging" by L. Xiao and S. Boyd
% "Fastest mixing Markov chain on a graph" by S. Boyd, P. Diaconis, and L. Xiao
% "Convex Optimization of Graph Laplacian Eigenvalues" by S. Boyd
%
% Almir Mutapcic 08/29/06

% maximum degree solution
[n,m] = size(A);

% max degrees of the nodes
Lunw = A*A';        % unweighted Laplacian matrix
degs = diag(Lunw);

% max degree weight allocation
max_deg = max(degs);
w = (1/max_deg)*ones(m,1);
