function [R1,R2] = makeSparseGradRemainderMatrices(numRows,numCols)

% Let D be the discrete gradient operator using symmetric BCs.
% D = G + R, where G uses periodic BCs and R is a sparse "leftover" matrix.
% R = [R1;R2].  This function computes R1 and R2.

numPix = numRows*numCols;

idxs = reshape(1:numPix,numRows,numCols);

rowIdxs = [];
colIdxs = [];
vals = [];

for col = 1:numCols
   
    rowIdx = idxs(numRows,col);
    colIdx = idxs(numRows,col);
    val = 1;
    
    rowIdxs = [rowIdxs,rowIdx];
    colIdxs = [colIdxs,colIdx];
    vals = [vals,val];
    
    rowIdx = idxs(numRows,col);
    colIdx = idxs(1,col);
    val = -1;
    
    rowIdxs = [rowIdxs,rowIdx];
    colIdxs = [colIdxs,colIdx];
    vals = [vals,val];        
    
end

R1 = sparse(rowIdxs,colIdxs,vals,numPix,numPix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rowIdxs = [];
colIdxs = [];
vals = [];

for row = 1:numRows
   
    rowIdx = idxs(row,numCols);
    colIdx = idxs(row,numCols);
    val = 1;
    
    rowIdxs = [rowIdxs,rowIdx];
    colIdxs = [colIdxs,colIdx];
    vals = [vals,val];
    
    rowIdx = idxs(row,numCols);
    colIdx = idxs(row,1);
    val = -1;
    
    rowIdxs = [rowIdxs,rowIdx];
    colIdxs = [colIdxs,colIdx];
    vals = [vals,val];        
    
end

R2 = sparse(rowIdxs,colIdxs,vals,numPix,numPix);

end