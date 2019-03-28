function out = prox2Norm(xHat,t)

nrm = norm(xHat(:));
if t <= nrm
    out = (1-t/nrm)*xHat;
else
    out = zeros(size(xHat));
end

end