function out = shearletTransform(x,qmf1,qmf2,scale,ndir)
% This uses SHEARLAB.
% In the Shearlab-1.1 directory, see the
% Shearlet_Transform_with_Separable_Generator directory.

[dd1 dd2] = sampled_DST(x,qmf1,qmf2,scale,ndir);

out = cat(3,dd1,dd2);

end