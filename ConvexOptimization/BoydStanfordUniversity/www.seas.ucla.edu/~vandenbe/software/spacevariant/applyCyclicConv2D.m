function out = applyCyclicConv2D(x,eigValArr)

% [numRows,numCols] = size(x);
% 
% a = zeros(numRows,numCols);
% a(1,1) = 1;
% 
% Ra = imfilter(a,mask,'circular');
% 
% RaHat = fft2(Ra);
% 
% eigValArr = conj(RaHat);

out = ifft2(eigValArr.*fft2(x));


end