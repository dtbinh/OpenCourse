% Floor planning example with an optimal trade-off curve.
% (a figure is generated)
%
% This is an example taken directly from the paper:
%
%   A Tutorial on Geometric Programming (see pages 24-25)
%   by Boyd, Kim, Vandenberghe, and Hassibi.
%
% Solves the problem of configuring and placing rectangles such
% that they do not overlap and that they minimize the area of the
% bounding box. This code solves the specific instances given
% in the GP tutorial. We have four rectangles with variable
% width w_i and height h_i. They need to satisfy area and aspect
% ration constraints. The GP is formulated as:
%
%   minimize   max(wa + wb, wc + wd)*(max(ha,hb) + max(hc,hd))
%       s.t.   wa*ha == area_a, wb*hb == area_b, ...
%              1/alpha_max <= ha/wa <= alpha_max, ...
%
% where variables are rectangle widths w's and heights h's.
%
% Almir Mutapcic 02/02/06

% constants
a = 0.2;
b = 0.5;
c = 1.5;
d = 0.5;

% GP variables
gpvar wa wb wc wd ha hb hc hd

% objective function is the area of the bounding box
obj = max(wa + wb, wc + wd)*(max(ha,hb) + max(hc,hd));

% constraints (now impose the non-changing constraints)
constr = [ ha*wa == a; hb*wb == b; hc*wc == c; hd*wd == d ];

% set the quiet flag (no solver reporting)
global QUIET; QUIET = 1;

% alpha is the changing parameter
N = 20;
alpha = linspace(1.01,4,N);

min_area = []; status_history = {};
for n = 1:N
  % add constraints that depend on the changing parameter alpha
  constr(5)  = 1/alpha(n) <= ha/wa; constr(6)  = ha/wa <= alpha(n);
  constr(7)  = 1/alpha(n) <= hb/wb; constr(8)  = hb/wb <= alpha(n);
  constr(9)  = 1/alpha(n) <= hc/wc; constr(10) = hc/wc <= alpha(n);
  constr(11) = 1/alpha(n) <= hd/wd; constr(12) = hd/wd <= alpha(n);

  [obj_value, solution, status] = gpsolve(obj, constr);
  min_area(n,1) = obj_value;
  status_history{end+1} = status;
end

% set the quiet flag (no solver reporting)
global QUIET; QUIET = 1;

plot(alpha,min_area);
xlabel('alpha'); ylabel('min area');
axis([1 4 2.5 4]);
