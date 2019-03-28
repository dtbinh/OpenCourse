% Exercise 4.31: Design of a cantilever beam (non-recursive formulation)
% (For a detailed explanation see section 4.5.4, pp. 163-165)
% Boyd & Vandenberghe "Convex Optimization"
% (a figure is generated)
%
% We have a segmented cantilever beam with N segments. Each segment
% has a unit length and variable width and height (rectangular profile).
% The goal is minimize the total volume of the beam, over all segment
% widths w_i and heights h_i, subject to constraints on aspect ratios,
% maximum allowable stress in the material, vertical deflection y, etc.
%
% The problem can be posed as a geometric program (posynomial form)
%     minimize   sum( w_i* h_i)
%         s.t.   w_min <= w_i <= w_max,       for all i = 1,...,N
%                h_min <= h_i <= h_max
%                S_min <= h_i/w_i <= S_max
%                6*i*F/(w_i*h_i^2) <= sigma_max
%                6*F/(E*w_i*h_i^3) == d_i
%                (2*i - 1)*d_i + v_(i+1) <= v_i
%                (i - 1/3)*d_i + v_(i+1) + y_(i+1) <= y_i
%                y_1 <= y_max
%
% with variables w_i, h_i, d_i, (i = 1,...,N) and v_i, y_i (i = 1,...,N+1).
% (Consult the book for other definitions and a recursive formulation of
% this problem.)
%
% Almir Mutapcic 01/25/06

% optimization variables
N = 8;
gpvar w(N) h(N) v(N+1) y(N+1);

% constants
wmin = .1; wmax = 100;
hmin = .1; hmax = 6;
Smin = 1/5; Smax = 5;
sigma_max = 1;
ymax = 10;
E = 1; F = 1;

% objective is the total volume of the beam
% obj = sum of (widths*heights*lengths) over each section
% (recall that the length of each segment is set to be 1)
obj = w'*h; 

% non-recursive formulation
d = 6*F*ones(N,1)./(E*ones(N,1).*w.*h.^3); 
constr_v = [];
constr_y = [];
for i = 1:N
  constr_v = [constr_v; (2*i-1)*d(i) + v(i+1) <= v(i)];
  constr_y = [constr_y; (i-1/3)*d(i) + v(i+1) + y(i+1) <= y(i)];
end

% constraint set
constr = [ ...
  wmin*ones(N,1) <= w; w <= wmax*ones(N,1);
  hmin*ones(N,1) <= h; h <= hmax*ones(N,1);
  Smin*ones(N,1) <= h./w; h./w <= Smax*ones(N,1);
  6*F*[1:N]'./(w.*(h.^2)) <= sigma_max*ones(N,1);
  constr_v; constr_y;
  y(1) <= ymax;
];

% solve GP and compute the optimal volume
[obj_value, solution, status] = gpsolve(obj, constr);
assign(solution);

% display results
disp('The optimal widths and heights are: ');
w, h
fprintf(1,'The optimal minimum volume of the beam is %3.4f\n', sum(w.*h))

% plot the 3D model of the optimal cantilever beam
close all;
plot_cbeam([h; w])
