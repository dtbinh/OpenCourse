function out = shearletAdjTransform(y,qmf1,qmf2,scale,ndir)

% This uses SHEARLAB.
% In the Shearlab-1.1 directory, see the
% Shearlet_Transform_with_Separable_Generator directory.

alpha = 2^(ndir+2)+2;

dd1 = y(:,:,1:alpha/2);
dd2 = y(:,:,alpha/2+1:end);

x1 = half_adj_sampled_DST1(dd1,qmf1,qmf2,scale,ndir,0);
x2 = half_adj_sampled_DST1(dd2,qmf1,qmf2,scale,ndir,1);

out = x1 + x2;

end