function out = projectOntoDualIsoNormUnitBall(y)

p1 = y(:,:,1);
p2 = y(:,:,2);

norms = sqrt(p1.^2 + p2.^2);

fnd = find(norms > 1);
p1(fnd) = p1(fnd)./norms(fnd);
p2(fnd) = p2(fnd)./norms(fnd);

out = cat(3,p1,p2);

end