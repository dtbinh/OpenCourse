function H = entropy(px)
% ENTROPY Computes the entropy of discrete r.v. with pmf px. 
%
%       ENTROPY computes the entropy of a discrete random variable
%       with probability mass function p_x.
%
% Almir Mutapcic 01/18/06
H = -sum( px.*log2(px) );
