function [xPlus_x1,costs] = deblur_mixedTV_decomp_drPrimal(b,kernels,Ks,params)
               
% This code corresponds to section 5.2 in the paper
% "Image deblurring by primal-dual operator splitting." 
% We solve minimize || Kx - b||_1 + gamma*|| Dx ||_{iso}
% subject to  0 <= x <= 1
% K = Kp + Ks. D = Dp + Ds.
% D uses symmetric BCs.

% In this code primal and dual step sizes are handled as follows:
% the operator A gets multiplied by beta, and g is replaced by
% gTilde, as in p. 7 of the paper.  
% With these modifications to the objective function,
% we then use Douglas-Rachford with the single step size t.

% After modifying K,D, and b, our problem is to 
% minimize (1/beta)||Kx - b||_1 + (gamma/beta)||Dx||_{iso}.

maxIter = params.maxIter;
t = params.t; % primal step size
beta = params.beta;
gamma = params.gamma;
showTrigger = params.showTrigger;
overRelax = params.overRelax; % called rho in paper.

evalIsoNorm = @(y) sum(sum(sqrt(y(:,:,1).^2 + y(:,:,2).^2)));

[numRows,numCols] = size(b);
numPix = numRows*numCols;

Ks = beta*Ks;
kernels = beta*kernels;
b = beta*b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nrMask,ncMask,R,S] = size(kernels);

nr = numRows/R;
nc = numCols/S;

eigValArrs_KpBlocks = zeros(nr,nc,R,S);
eigValArrs_KpTransBlocks = zeros(nr,nc,R,S);
for r = 1:R
    for s = 1:S
        kernel = kernels(:,:,r,s);
        arr = eigValArrForCyclicConvOp(kernel,nr,nc);
        eigValArrs_KpBlocks(:,:,r,s) = arr;
        eigValArrs_KpTransBlocks(:,:,r,s) = conj(arr);
    end
end

applyKp = @(x) applyBlockCyclicConv2D(x,eigValArrs_KpBlocks);
applyKpTrans = @(x) applyBlockCyclicConv2D(x,eigValArrs_KpTransBlocks);

applyKs = @(x) reshape(Ks*x(:),numRows,numCols);
applyKsTrans = @(x) reshape(Ks'*x(:),numRows,numCols);

applyK = @(x) applyKp(x) + applyKs(x);
applyKTrans = @(x) applyKpTrans(x) + applyKsTrans(x);

Ds = formSparseMatrix_GradMinusBlockGrad_periodic(numRows,numCols,R,S); % Each block uses periodic BCs.
Ds1 = Ds(1:numPix,:);
Ds2 = Ds(numPix+1:end,:);

[rem1,rem2] = makeSparseGradRemainderMatrices(numRows,numCols);
Ds1 = Ds1 + rem1;
Ds2 = Ds2 + rem2;

Ds1 = beta*Ds1;
Ds2 = beta*Ds2;
applyDs1 = @(x) reshape(Ds1*x(:),numRows,numCols);
applyDs1Trans = @(x) reshape(Ds1'*x(:),numRows,numCols);
applyDs2 = @(x) reshape(Ds2*x(:),numRows,numCols);
applyDs2Trans = @(x) reshape(Ds2'*x(:),numRows,numCols);

applyDs = @(x) cat(3,applyDs1(x),applyDs2(x));
applyDsTrans = @(y) applyDs1Trans(y(:,:,1)) + applyDs2Trans(y(:,:,2));

eigValArr_D1Block = eigValArrForCyclicConvOp(beta*[-1 1]',nr,nc);
eigValArr_D2Block = eigValArrForCyclicConvOp(beta*[-1 1],nr,nc);

eigValArr_D1BlockTrans = conj(eigValArr_D1Block);
eigValArr_D2BlockTrans = conj(eigValArr_D2Block);

applyD1Block = @(x) applyCyclicConv2D(x,eigValArr_D1Block);
applyD1BlockTrans = @(x) applyCyclicConv2D(x,eigValArr_D1BlockTrans);

applyD2Block = @(x) applyCyclicConv2D(x,eigValArr_D2Block);
applyD2BlockTrans = @(x) applyCyclicConv2D(x,eigValArr_D2BlockTrans);

% applyD1Block = @(x) beta*computeGradx_neumannBCs(x);
% applyD2Block = @(x) beta*computeGrady_neumannBCs(x);

applyDp1 = @(x) applyOpToEachImgBlock(x,applyD1Block,R,S);
applyDp1Trans = @(x) applyOpToEachImgBlock(x,applyD1BlockTrans,R,S);
applyDp2 = @(x) applyOpToEachImgBlock(x,applyD2Block,R,S);
applyDp2Trans = @(x) applyOpToEachImgBlock(x,applyD2BlockTrans,R,S);

applyDp = @(x) cat(3,applyDp1(x),applyDp2(x));
applyDpTrans = @(y) applyDp1Trans(y(:,:,1)) + applyDp2Trans(y(:,:,2));

applyD = @(x) applyDp(x) + applyDs(x);
applyDTrans = @(y) applyDpTrans(y) + applyDsTrans(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

applyB = @(x) cat(3,applyKp(x),applyDp(x));
applyBTrans = @(y) applyKpTrans(y(:,:,1)) + applyDpTrans(y(:,:,2:end));

applyC = @(x) cat(3,applyKs(x),applyDs(x));
applyCTrans = @(y) applyKsTrans(y(:,:,1)) + applyDsTrans(y(:,:,2:end));

applyA = @(x) cat(3,applyK(x),applyD(x));
applyATrans = @(y) applyKTrans(y(:,:,1)) + applyDTrans(y(:,:,2:end));

applyMtrx = @(x) x + applyKpTrans(applyKp(x)) + applyDpTrans(applyDp(x));
% eigValsMtrx = eigValArr_KpTrans.*eigValArr_Kp + eigValArr_Dp1Trans.*eigValArr_Dp1 ...
%                 + eigValArr_Dp2Trans.*eigValArr_Dp2 + ones(numRows,numCols);
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eigValsMtrxBlocks = zeros(nr,nc,R,S);
for r = 1:R
    for s = 1:S
        
        eigValArr_KpBlock = eigValArrs_KpBlocks(:,:,r,s);
        eigValArr_KpTransBlock = eigValArrs_KpTransBlocks(:,:,r,s);
        
        eigValsMtrxBlocks(:,:,r,s) = eigValArr_KpTransBlock.*eigValArr_KpBlock + ...
                                eigValArr_D1BlockTrans.*eigValArr_D1Block + eigValArr_D2BlockTrans.*eigValArr_D2Block  ...
                                + ones(nr,nc);         
                            
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
sparseMtrx = speye(numPix) + Ks'*Ks + Ds1'*Ds1 + Ds2'*Ds2;

[U,statusInfo,perm] = chol(sparseMtrx,'vector');
UTrans = U';
permInv(perm) = 1:numPix;

p_x1 = b;
p_x2 = b;
p_y1 = applyB(b);
p_y2 = applyC(b);

costs = [];
figure('Name','xPlus_x inside DR')

for n = 1:maxIter 

    % Evaluate resolvent at p_x1,p_x2,p_y1,p_y2.
    xPlus_x1 = min((p_x1 + p_x2)/2,1);
    xPlus_x1 = max(xPlus_x1,0);
    xPlus_x2 = xPlus_x1;
    
    term = p_y1 + p_y2; % z is prox op. of g evaluated at term.
    termData = term(:,:,1);
    termReg = term(:,:,2:end);
    
    zData = prox1Norm(termData - b,2*t/beta) + b;
    zReg = proxIsoNorm(2*t*gamma/beta,termReg);
    z = cat(3,zData,zReg);    
    xPlus_y1 = (z + p_y1 - p_y2)/2;
    xPlus_y2 = (z + p_y2 - p_y1)/2;
    
    x1Hat = 2*xPlus_x1 - p_x1;
    x2Hat = 2*xPlus_x2 - p_x2;
    y1Hat = 2*xPlus_y1 - p_y1;
    y2Hat = 2*xPlus_y2 - p_y2;
    
    rhs = x1Hat + applyBTrans(y1Hat); %%
        
    yPlus_x1 = zeros(numRows,numCols);    
    for r = 1:R
        for s = 1:S
            
            rowStart = nr*(r-1) + 1;
            rowStop = nr*(r-1) + nr;
            colStart = nc*(s-1) + 1;
            colStop = nc*(s-1) + nc;   
            
            rhsBlock = rhs(rowStart:rowStop,colStart:colStop);
            yPlus_x1(rowStart:rowStop,colStart:colStop) = ifft2(fft2(rhsBlock)./eigValsMtrxBlocks(:,:,r,s));
            
        end
    end
    
%     yPlus_x1 = ifft2(fft2(rhs)./eigValsMtrx);
    yPlus_y1 = applyB(yPlus_x1);
    
    rhs2 = x2Hat + applyCTrans(y2Hat);
    rhs2 = rhs2(:); 
    temp = UTrans\rhs2(perm);    
    w = U\temp;    
    yPlus_x2 = w(permInv);
    yPlus_x2 = reshape(yPlus_x2,[numRows,numCols]);   
    yPlus_y2 = applyC(yPlus_x2);
    
%     check = sparseMtrx*yPlus_x2(:) - rhs2(:);
%     check = max(abs(check(:)));
    
    pPlus_x1 = p_x1 + overRelax*(yPlus_x1 - xPlus_x1);
    pPlus_x2 = p_x2 + overRelax*(yPlus_x2 - xPlus_x2);
    pPlus_y1 = p_y1 + overRelax*(yPlus_y1 - xPlus_y1);
    pPlus_y2 = p_y2 + overRelax*(yPlus_y2 - xPlus_y2);
    
    p_x1 = pPlus_x1;
    p_x2 = pPlus_x2;
    p_y1 = pPlus_y1;
    p_y2 = pPlus_y2;
               
    Dx = applyD(xPlus_x1);
    KxMinusb = applyK(xPlus_x1) - b;
    primalCost = (1/beta)*sum(abs(KxMinusb(:))) + (gamma/beta)*evalIsoNorm(Dx);
    
    costs = [costs,primalCost];         
                    
    if mod(n,showTrigger) == 0         
        imshow(xPlus_x1,[])    
        disp(['drPrimal iteration is: ',num2str(n)])
        disp(['primal cost is: ',num2str(primalCost)])         
        keyboard        
    end                 
                
end


end