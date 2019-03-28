function mtrx = formSparseMatrix_GradMinusBlockGrad_periodic(numRows,numCols,R,S)
% This code assumes each block uses PERIODIC boundary conditions.
% To understand this code, there's an output for each pixel in the image.
% (Actually two outputs, because gradient has two components.)
% Think about which outputs need to be modified.
% This code creates a matrix that modifies the appropriate outputs.

numPix = numRows*numCols;

nr = numRows/R;
nc = numCols/S;

numNonZero = (R-1)*numCols*2 + (S-1)*numRows*2;

mtrx = sparse([],[],[],2*numRows*numCols,numRows*numCols,numNonZero);

for r = 1:(R-1)
        
    i = nr*(r-1) + nr;                
        
    for j = 1:numCols
       
        ijIdx = (j-1)*numRows + i;
        ip1jIdx = (j-1)*numRows + i + 1;
        
        rowIdx = ijIdx;
%         colIdx = ijIdx;
%         mtrx(rowIdx,colIdx) = -1;
        
        colIdx = ip1jIdx;
        mtrx(rowIdx,colIdx) = 1;
        
%         ip1_mod = nr*(r-1) + 1; % = i + 1 - nr
        ip1_mod = i + 1 - nr;
        colIdx = (j-1)*numRows + ip1_mod;
        mtrx(rowIdx,colIdx) = -1;
                        
    end
       
end

i = numRows;
for j = 1:numCols
    
    ijIdx = (j-1)*numRows + i;        
    rowIdx = ijIdx;
    
    ip1jIdx = (j-1)*numRows + 1;
    colIdx = ip1jIdx;
    mtrx(rowIdx,colIdx) = 1;
    
    ip1_mod = numRows + 1 - nr;
    colIdx = (j-1)*numRows + ip1_mod;
    mtrx(rowIdx,colIdx) = -1;
    
end

for s = 1:(S-1)
    
    j = nc*(s-1) + nc;         
    
    for i = 1:numRows
        
        ijIdx = (j-1)*numRows + i;
        ijp1Idx = (j+1-1)*numRows + i;
        
        rowIdx = numPix + ijIdx;
%         colIdx = ijIdx;        
%         mtrx(rowIdx,colIdx) = -1;
        
        colIdx = ijp1Idx;
        mtrx(rowIdx,colIdx) = 1;
        
        jp1_mod = j + 1 - nc;
        colIdx = (jp1_mod - 1)*numRows + i;
        mtrx(rowIdx,colIdx) = -1;
        
    end
    
end

j = numCols;
for i = 1:numRows
    
    ijIdx = (j-1)*numRows + i;   
    rowIdx = numPix + ijIdx;
    
    jp1 = 1;
    ijp1Idx = (jp1 - 1)*numRows + i; % = i
    colIdx = ijp1Idx;
    mtrx(rowIdx,colIdx) = 1;
    
    jp1_mod = numCols + 1 - nc;
    colIdx = (jp1_mod - 1)*numRows + i;
    mtrx(rowIdx,colIdx) = -1;
    
end


end