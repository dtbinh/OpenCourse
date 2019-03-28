function eigValArr = eigValArrForCyclicConvOp(mask,numRows,numCols)

a = zeros(numRows,numCols);
a(1,1) = 1;
Ra = imfilter(a,mask,'circular');
RaHat = fft2(Ra);
eigValArr = RaHat;
% eigValArr = conj(RaHat);


end