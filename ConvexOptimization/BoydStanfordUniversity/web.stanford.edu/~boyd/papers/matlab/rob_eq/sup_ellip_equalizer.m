function cvx_optval = sup_ellip_equalizer(h,w,G_des,G_nom,p,q)
% SUP_ELLIP_EQUALIZER robust equalizer over an ellipsoid.
%
% Gives an epigraph formulation of the following convex function g(h):
%
%   g(h) = sup_{ G \in E } |G*w*h - G_des|
%
% where ellipsoid E = { G_nom + p*u1 + q*u2 | ||u|| <= 1 }
% and p, q are complex numbers, and u = (u1, u2) is a real vector.

cvx_begin
  epigraph variable t;
  variable lambda;
    lambda >= 0;
    [ lambda     0       0        real(p*w*h)   imag(p*w*h);
      0       lambda     0        real(q*w*h)   imag(q*w*h);
      0          0    t-lambda    real(G_nom*w*h-G_des)   imag(G_nom*w*h-G_des);
      real(p*w*h)  real(q*w*h)  real(G_nom*w*h-G_des)   t   0;
      imag(p*w*h)  imag(q*w*h)  imag(G_nom*w*h-G_des)   0   t;
    ] == hermitian_semidefinite(5);
cvx_end
