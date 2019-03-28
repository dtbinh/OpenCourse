function [w] = mh(A)
% MH	Returns the Metropolis-Hastings heuristic edge weights.
%
% [w] = MH(A) gives a vector of the Metropolis-Hastings edge weights
% for a graph described by the incidence matrix A (n x m).
% Here n is the number of nodes and m is the number of edges in the graph;
% each column of A has exactly one +1 and one -1.
%
% M.-H. weight on an edge is one over the maximum of the degrees
% of the adjacent nodes.
%
% For more details see the references:
% "Fast linear iterations for distributed averaging" by L. Xiao and S. Boyd
% "Fastest mixing Markov chain on a graph" by S. Boyd, P. Diaconis, and L. Xiao
% "Convex Optimization of Graph Laplacian Eigenvalues" by S. Boyd
%
% Almir Mutapcic 08/29/06

% degrees of the nodes
Lunw = A*A';          % unweighted Laplacian matrix
degs = diag(Lunw);

% Metropolis-Hastings weights
mh_degs = abs(A)'*diag(degs);
w = 1./max(mh_degs,[],2);
