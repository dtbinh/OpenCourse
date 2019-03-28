function [uk,costs] = deblur_mixedTV_decomp_admm(b,kernels,Ks,params)
               
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
beta = params.beta; % for scaling the problem.
gamma = params.gamma; 
showTrigger = params.showTrigger;
alpha = params.overRelax; 

evalIsoNorm = @(y) sum(sum(sqrt(y(:,:,1).^2 + y(:,:,2).^2)));

[numRows,numCols] = size(b);
numPix = numRows*numCols;

Ks = beta*Ks;
kernels = beta*kernels;
b = beta*b;

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

x1k = b;
x2k = b;
uk = b;
y1k = applyB(x1k);
y2k = applyC(x2k);

w1k = zeros(size(uk));
w2k = zeros(size(uk));
z1k = zeros(size(y1k));
z2k = zeros(size(y2k));

costs = [];

figure('Name','uk inside ADMM')

for k = 1:maxIter
    
    % Minimize aug. Lagrangian w.r.t. (x1,x2).
    rhs1 = uk - w1k/t + applyBTrans(y1k - z1k/t);
%     x1kp1 = ifft2(fft2(rhs1)./eigValsMtrx);

    x1kp1 = zeros(numRows,numCols);    
    for r = 1:R
        for s = 1:S
            
            rowStart = nr*(r-1) + 1;
            rowStop = nr*(r-1) + nr;
            colStart = nc*(s-1) + 1;
            colStop = nc*(s-1) + nc;   
            
            rhsBlock = rhs1(rowStart:rowStop,colStart:colStop);
            x1kp1(rowStart:rowStop,colStart:colStop) = ifft2(fft2(rhsBlock)./eigValsMtrxBlocks(:,:,r,s));
            
        end
    end

%     check = applyMtrx(x1kp1) - rhs1;
%     check = max(abs(check(:)));
    
    rhs2 = uk - w2k/t + applyCTrans(y2k - z2k/t);    
    rhs2 = rhs2(:); 
    temp = UTrans\rhs2(perm);    
    q = U\temp;    
    x2kp1 = q(permInv);
    x2kp1 = reshape(x2kp1,[numRows,numCols]);   
%     check = sparseMtrx*x2kp1(:) - rhs2;
%     check = max(abs(check(:)));
    
    % Minimize aug. Lagrangian w.r.t. (u,y1,y2);
    
    term = .5*(x1kp1 + w1k/t + x2kp1 + w2k/t);
    ukp1 = min(term,1);
    ukp1 = max(ukp1,0);
    
    y1Hat = applyB(x1kp1) + z1k/t;
    y2Hat = applyC(x2kp1) + z2k/t;
    
    term = y1Hat + y2Hat;
    termData = term(:,:,1);
    termReg = term(:,:,2:end);
    
    tmpData = prox1Norm(termData - b,2/(t*beta)) + b;
    tmpReg = proxIsoNorm(2*gamma/(t*beta),termReg);
    tmp = cat(3,tmpData,tmpReg);    
    
    y1kp1 = (tmp + y1Hat - y2Hat)/2;
    y2kp1 = (tmp + y2Hat - y1Hat)/2;
    
%     w1kp1 = w1k + t*(x1kp1 - ukp1);
%     w2kp1 = w2k + t*(x2kp1 - ukp1);
%     z1kp1 = z1k + t*(applyB(x1kp1) - y1kp1);
%     z2kp1 = z2k + t*(applyC(x2kp1) - y2kp1);    
    
    w1kp1 = w1k + t*(alpha*x1kp1 + (1-alpha)*uk - ukp1);
    w2kp1 = w2k + t*(alpha*x2kp1 + (1-alpha)*uk - ukp1);
    z1kp1 = z1k + t*(alpha*applyB(x1kp1) + (1-alpha)*y1k - y1kp1);
    z2kp1 = z2k + t*(alpha*applyC(x2kp1) + (1-alpha)*y2k - y2kp1);     
    
    x1k = x1kp1;
    x2k = x2kp1;
    uk = ukp1;
    y1k = y1kp1;
    y2k = y2kp1;
    
    w1k = w1kp1;
    w2k = w2kp1;
    z1k = z1kp1;
    z2k = z2kp1;
    
    Dx = applyD(uk);
    KxMinusb = applyK(uk) - b;
    primalCost = (1/beta)*sum(abs(KxMinusb(:))) + (gamma/beta)*evalIsoNorm(Dx);
    
    costs = [costs,primalCost];         
                    
    if mod(k,showTrigger) == 0         
        imshow(uk,[])    
        disp(['admm iteration is: ',num2str(k)])
        disp(['primal cost is: ',num2str(primalCost)])         
        keyboard        
    end         
    
    
end

end















