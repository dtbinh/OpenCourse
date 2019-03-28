function [y, yp, ypp] = bsplines(x)

% [y, yp, ypp] = bsplines(x)
%
% Computes the values and the first and second derivatives of the 13 cubic 
% B-splines with breakpoints 0, 1, 2, ..., 10.   
% x must be in the interval [0, 10].
%
% y:    a column vector with the values of the 13 cubic B-splines evaluated
%       at x.
% yp:   a column vector with the first derivatives of the 13 cubic B-splines
%       B-splines evaluated at x.
% ypp:  a column vector with the second derivatives of the 13 cubic 
%       B-splines evaluated at x.

if x < 0 || x > 10
   error('x must be in [0, 10]');
end;

% knots
a = [ 0, 0, 0, 0:10, 10, 10, 10 ];

% 0 <= x < 10 and a(i) <= x < a(i+1) or x = 10 and i is 14.
i = floor(x) + 4;

% y1, y2, y are vectors of linear, quadratic, and cubic B-splines at x.
if i == 14, 
    y1 = [ zeros(10,1); 1 ];
    y2 = [ zeros(11,1); 1 ];
    y = [ zeros(12,1); 1 ];
else
    w1 = (x - a(i)) / (a(i+1) - a(i));
    G1 = [ 1 - w1; w1 ];
    w2 = [ (x - a(i-1)) / (a(i+1) - a(i-1) ) ; 
           (x - a(i)) / (a(i+2) - a(i)) ];
    G2 = [ diag(1-w2); 0, 0] + [ 0, 0; diag(w2) ];
    w3 = [ (x - a(i-2)) / (a(i+1) - a(i-2)) ;
           (x - a(i-1)) / (a(i+2) - a(i-1)) ;
           (x - a(i)) / (a(i+3) - a(i)) ];
    G3 = [ diag(1-w3); 0, 0, 0 ] + [ 0, 0, 0; diag(w3) ];
    y = zeros(13, 1);
    y(i-3:i) = G3*G2*G1; 
    y2 = zeros(12, 1);
    y2(i-3:i-1) = G2*G1;
    y1 = zeros(11, 1);
    y1(i-3:i-2) = G1;
end;

% Width of support of linear, quadratic, and cubic B-splines.
d1 = a(5:15)' - a(3:13)';  
d2 = a(5:16)' - a(2:13)';  
d3 = a(5:17)' - a(1:13)';  

% Derivatives of quadratic B-splines.
y2p = zeros(12,1);
y2p(1) = -2 * y1(1) / d1(1);
y2p(2:11) = 2 * ( y1(1:10) ./ d1(1:10) - y1(2:11) ./ d1(2:11) );
y2p(12) = 2 * y1(11) / d1(11);

% First derivatives of cubic B-splines.
yp = zeros(13,1);
yp(1) = -3 * y2(1) / d2(1);
yp(2:12) = 3 * ( y2(1:11) ./ d2(1:11) - y2(2:12) ./ d2(2:12) );
yp(13) = 3 * y2(12) / d2(12);

% Second derivatives of cubic B-splines.
ypp = zeros(13,1);
ypp(1) = -3 * y2p(1) / d2(1);
ypp(2:12) = 3 * ( y2p(1:11) ./ d2(1:11) - y2p(2:12) ./ d2(2:12) );
ypp(13) = 3 * y2p(12) / d2(12);
