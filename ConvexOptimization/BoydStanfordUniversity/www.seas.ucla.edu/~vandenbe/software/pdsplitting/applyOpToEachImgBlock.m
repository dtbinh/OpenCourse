function out = applyOpToEachImgBlock(x,applyOp,R,S)

[numRows,numCols] = size(x);
nr = numRows/R;
nc = numCols/S;

out = zeros(numRows,numCols);

for r = 1:R
    for s = 1:S

        rowStart = nr*(r-1) + 1;
        rowStop = nr*(r-1) + nr;
        colStart = nc*(s-1) + 1;
        colStop = nc*(s-1) + nc;
        
        block = x(rowStart:rowStop,colStart:colStop);
        out(rowStart:rowStop,colStart:colStop) = applyOp(block);
        
    end
end        

end