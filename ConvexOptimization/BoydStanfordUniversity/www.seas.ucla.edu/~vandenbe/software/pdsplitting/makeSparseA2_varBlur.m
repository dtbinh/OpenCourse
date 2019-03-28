function A2 = makeSparseA2_varBlur(kernels,numRows,numCols)

% A1 blurs each block independently, using periodic BCs on each block.
% A blurs entire image, using different kernel for each block, and periodic
% BCs at the boundary of the image.
% A = A1 + A2.

% All kernels have the same size.

[nrMask,ncMask,R,S] = size(kernels); % This code assumes nrMask and ncMask are odd.

numPix = numRows*numCols;

halfr = (nrMask - 1)/2;
halfc = (ncMask - 1)/2;

nr = numRows/R;
nc = numCols/S;

innerBlockArea = (nr - 2*halfr)*(nc - 2*halfc);
blockBorderArea = nr*nc - innerBlockArea;
numNonZero = blockBorderArea*R*S*2*nrMask*ncMask;

rowArr = zeros(numNonZero,1);
colArr = zeros(numNonZero,1);
valArr = zeros(numNonZero,1);

count = 1;

% border = zeros(numRows,numCols);
% for r = 1:R
%     for s = 1:S
%         
%         rowStart = nr*(r-1) + 1;
%         rowStop = nr*(r-1) + nr;
%         colStart = nc*(s-1) + 1;
%         colStop = nc*(s-1) + nc;
%         
%         border(rowStart + halfr:(rowStop - halfr),...
%                colStart + halfc:(colStop - halfc)) = 1;        
%         
%         
%     end
% end
% border = 1 - border;
% figure('Name','border')
% colormap(gray)
% imagesc(border)

for r = 1:R
    for s = 1:S
        
        mask = kernels(:,:,r,s);

        rowStart = nr*(r-1) + 1;
        rowStop = nr*(r-1) + nr;
        colStart = nc*(s-1) + 1;
        colStop = nc*(s-1) + nc;
        
        border = zeros(nr,nc);
        border(1 + halfr:(nr - halfr),...
               1 + halfc:(nc - halfc)) = 1; 
        
        border = 1 - border;
        arr = zeros(numRows,numCols);
        arr(rowStart:rowStop,colStart:colStop) = border;
        
        [rows,cols] = find(arr);
        
        for idx = 1:length(rows)
            
            iBase = rows(idx);
            jBase = cols(idx); 
            
%             mask = kernels(iBase,jBase,:);
%             mask = reshape(mask,[nrMask,ncMask]);
            
            rowIdx = numRows*(jBase - 1) + iBase;            
            
            for dj = -halfc:halfc
                for di = -halfr:halfr 
                    
                    i = iBase + di;
                    j = jBase + dj;                    
                    
                    if (i < rowStart) | (i > rowStop) | (j < colStart) | (j > colStop)  
                        
                        iModOuter = mod(i-1,numRows) + 1;
                        jModOuter = mod(j-1,numCols) + 1;
                        
                        iMod = mod(i-1,nr) + rowStart;
                        jMod = mod(j-1,nc) + colStart;
                        
                        coeff = mask(di + halfr + 1, dj + halfc + 1);
                        
                        colIdx = numRows*(jMod - 1) + iMod;
                        rowArr(count) = rowIdx;
                        colArr(count) = colIdx;
                        valArr(count) = -coeff;
                        count = count + 1;
                        
                        colIdx = numRows*(jModOuter - 1) + iModOuter;
                        rowArr(count) = rowIdx;
                        colArr(count) = colIdx;
                        valArr(count) = coeff;
                        count = count + 1;                                                                            

                    end                    
                                                            
                end
            end
                                    
        end

    end
end

rowArr = rowArr(1:count-1);
colArr = colArr(1:count-1);
valArr = valArr(1:count-1);
numNonZero = length(valArr);

A2 = sparse(rowArr,colArr,valArr,numPix,numPix);


end