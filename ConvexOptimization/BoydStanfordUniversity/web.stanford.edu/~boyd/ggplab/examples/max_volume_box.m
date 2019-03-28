% Box volume maximization example.
% (a figure is generated if the tradeoff flag is turned on)
%
% This is an example taken directly from the paper:
%
%   A Tutorial on Geometric Programming (see pages 8 and 13)
%   by Boyd, Kim, Vandenberghe, and Hassibi.
%
% Maximizes volume of a box-shaped structure which has constraints
% on its total wall area, its total floor area, and which has lower
% and upper bounds on the aspect ratios. This leads to a GP:
%
%   maximize   h*w*d
%       s.t.   2(h*w + h*d) <= Awall, w*d <= Afloor
%              alpha <= h/w <= beta
%              gamma <= d/w <= delta
%
% where variables are the box height h, width w, and depth d.
%
% Almir Mutapcic 02/01/2006
clear all; close all;
PLOT_TRADEOFF = 1; % to disable tradeoff plot, set PLOT_TRADEOFF = 0

% problem constants
Awall  = 10000; Afloor = 1000;
alpha = 0.5; beta = 2; gamma = 0.5; delta = 2;

% GP variables
gpvar h w d

% objective function is the box volume
volume = h*w*d;

% set of constraints expressed as an array
constr = [ 2*(h*w + h*d) <= Awall; w*d <= Afloor;
           alpha <= h/w; h/w <= beta;
           gamma <= d/w; d/w <= delta;];

% solve the given GP problem
[max_volume solution status] = gpsolve(volume, constr, 'max');
assign(solution);

% display results
fprintf(1,'\nMaximum volume is %2.2f.\n', max_volume)
fprintf(1,['Optimal solution is height h = %3.4f, width w = %3.4f, '...
           'depth d = %3.4f.\n\n'], h, w, d);

%*********************************************************************
% tradeoff curve code
%*********************************************************************
if( PLOT_TRADEOFF )

% set the quiet flag (no solver reporting)
global QUIET; QUIET = 1;

% reintroduce the optimization variables
gpvar h w d

% varying parameters for an optimal trade-off curve
N = 10;
Aflr = logspace(1,3,N);
Awall = [100 1000 10000];
opt_volumes = zeros(length(Awall),N);

% setup various GP problems with varying parameters
for k = 1:length(Awall)
  for n = 1:N
    % change constraints with varying parameters
    constr(1) = 2*h*w + 2*h*d <= Awall(k);
    constr(2) = w*d <= Aflr(n);

    % solve the GP problem and compute the optimal volume
    [max_volume solution status] = gpsolve(volume, constr, 'max');
    opt_volumes(k,n) = max_volume;
  end
end

% enable solver reporting again
global QUIET; QUIET = 0;

% plot the tradeoff curve
loglog(Aflr,opt_volumes(1,:), Aflr,opt_volumes(2,:), Aflr,opt_volumes(3,:));
xlabel('Aflr'); ylabel('V');

end
% end tradeoff curve code
