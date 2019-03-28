function [aalpha] = aalpha(T)
% A helper function for simplified spot rate curve extraction.
% Computes value of the a_alpha exponent.
%
global beta
aalpha = (exp(beta*T) - 1)/beta^2 - T/beta;
aalpha = -aalpha;
