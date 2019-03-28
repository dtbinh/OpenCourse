function [cvx_optval, s_equal] = opt_equal_speed(G,Tother,Tamb,Tmax,smin,smax)
%
% suboptimal solution (find optimal constant speed for all processors)
%
[m n] = size(G);
cvx_begin
  variable s_equal    % processor speed variables
  maximize n*s_equal  % maximize throughput
  subject to
    % maximum temperature constraint
    G*ones(n,1)*pow_pos(s_equal,3) + Tother + Tamb <= Tmax;
    % processor speed bounds
    smin <= s_equal; s_equal <= smax;
cvx_end
