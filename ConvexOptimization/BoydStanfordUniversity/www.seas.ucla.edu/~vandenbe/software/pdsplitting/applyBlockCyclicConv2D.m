function out = applyBlockCyclicConv2D(x,eigValArrs)

[nr,nc,R,S] = size(eigValArrs);
out = zeros(size(x));
for r = 1:R
    for s = 1:S
        
        rowStart = nr*(r-1) + 1;
        rowStop = nr*(r-1) + nr;
        colStart = nc*(s-1) + 1;
        colStop = nc*(s-1) + nc;    
        
        arr = eigValArrs(:,:,r,s);
        out(rowStart:rowStop,colStart:colStop) = ifft2(arr.*fft2(x(rowStart:rowStop,colStart:colStop)));
        
    end
end

end