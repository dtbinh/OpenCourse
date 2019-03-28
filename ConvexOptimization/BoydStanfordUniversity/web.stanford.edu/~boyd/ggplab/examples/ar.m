function [ar] = ar(T)
% A helper function for simplified spot rate curve extraction.
% Computes value of the a_r exponent.
%
global beta
ar = (exp(beta*T) - 1)/beta;
ar = -ar;
