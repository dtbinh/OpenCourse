function A2 = makeSparseA2_replicateBCs_varBlur(numRows,numCols,kernels)
% I'll assume nrMask and ncMask are odd.

% [nrMask,ncMask] = size(mask);
[nrMask,ncMask,R,S] = size(kernels);

nr = numRows/R;
nc = numCols/S;

halfr = (nrMask - 1)/2;
halfc = (ncMask - 1)/2;

border = zeros(numRows,numCols);
border(nrMask:(numRows - nrMask + 1),ncMask:(numCols - ncMask+1)) = 1;
border = 1 - border;
[rows,cols] = find(border);

numNonZero = length(rows)*nrMask*ncMask + max(nrMask,ncMask);

A2 = sparse([],[],[],numRows*numCols,numRows*numCols,numNonZero);
        
for idx = 1:length(rows)
    
    iBase = rows(idx);
    jBase = cols(idx);
    
    r = ceil(iBase/nr);
    s = ceil(jBase/nc);
    mask = kernels(:,:,r,s);
            
    rowIdx = numRows*(jBase - 1) + iBase;
    
        for dj = -halfc:halfc
            for di = -halfr:halfr
                
                i = iBase + di;
                j = jBase + dj;
                
                if (i < 1) | (i > numRows) | (j < 1) | (j > numCols)                                       
                    
                    iMod = mod(i-1,numRows) + 1;
                    jMod = mod(j-1,numCols) + 1;                                        
                    
                    coeff = mask(di + halfr + 1, dj + halfc + 1);
%                     out(iBase,jBase) = out(iBase,jBase) - coeff*x(iMod,jMod);                    
                    
                    colIdx = numRows*(jMod - 1) + iMod;
                    A2(rowIdx,colIdx) = A2(rowIdx,colIdx) - coeff;
                    
                    if i < 1                        
                        closestRow = 1;                        
                    else
                        closestRow = min(i,numRows);
                    end
                    
                    if j < 1
                        closestCol = 1;
                    else
                        closestCol = min(j,numCols);
                    end                     
                    
%                     out(iBase,jBase) = out(iBase,jBase) + coeff*x(closestRow,closestCol);
                    
                    colIdx = numRows*(closestCol - 1) + closestRow;
                    A2(rowIdx,colIdx) = A2(rowIdx,colIdx) + coeff;                    
                    
                end
                
            end
        end 
                        
end
        

end