function [b] = b(T)
% A helper function for simplified spot rate curve extraction.
% Computes value of the b exponent.
%
global beta
b = (exp(beta*T) - 1)/beta^2 - T/beta;
b = -b;
