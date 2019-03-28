function out = applyNagyBlur(x,eigValArrs,U)

out = zeros(size(x));

numKernels = size(eigValArrs,3);

for k = 1:numKernels
        
    out = out + U(:,:,k).*applyCyclicConv2D(x,eigValArrs(:,:,k));

end


end