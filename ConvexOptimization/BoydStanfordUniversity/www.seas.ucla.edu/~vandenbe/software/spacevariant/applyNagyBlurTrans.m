function out = applyNagyBlurTrans(x,eigValArrs_KpTrans,U)

out = zeros(size(x));

numKernels = size(eigValArrs_KpTrans,3);

for k = 1:numKernels
    
    out = out + applyCyclicConv2D(U(:,:,k).*x,eigValArrs_KpTrans(:,:,k));    
    
end


end