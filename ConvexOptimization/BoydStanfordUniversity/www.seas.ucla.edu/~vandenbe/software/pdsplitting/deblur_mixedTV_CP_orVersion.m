function [xn,costs] = deblur_mixedTV_CP_orVersion(b,mask,S,params)
% orVersion stands for "overrelaxed version" of Chambolle-Pock.

% We solve minimize ||Ax - b|| + gamma*|| Gx ||
% subject to 0 <= x_{ij} <= 1 for all i,j.
% A convolves with mask using replicate boundary conditions.
% G is a discrete gradient operator using symmetric boundary conditions.
% A = P + S , G = D + R.

maxIter = params.maxIter;
showTrigger = params.showTrigger;
gamma = params.gamma;
theta = params.theta;
overRelax = params.overRelax;

evalIsoNorm = @(y) sum(sum(sqrt(y(:,:,1).^2 + y(:,:,2).^2)));

[numRows,numCols] = size(b);
numPix = numRows*numCols;
%%%%%%%%%%%% Set up operators.
adjMask = fliplr(mask);
adjMask = flipud(adjMask);

eigValArr_P = eigValArrForCyclicConvOp(mask,numRows,numCols);
eigValArr_PTrans = conj(eigValArr_P);

applyP = @(x) applyCyclicConv2D(x,eigValArr_P);
applyPTrans = @(x) applyCyclicConv2D(x,eigValArr_PTrans);

applyS = @(x) reshape(S*x(:),numRows,numCols);
applySTrans = @(x) reshape(S'*x(:),numRows,numCols);

applyA = @(x) applyP(x) + applyS(x);
applyATrans = @(x) applyPTrans(x) + applySTrans(x);

%%%%

eigValArr_D1 = eigValArrForCyclicConvOp([-1 1]',numRows,numCols);
eigValArr_D2 = eigValArrForCyclicConvOp([-1 1],numRows,numCols);

eigValArr_D1Trans = conj(eigValArr_D1);
eigValArr_D2Trans = conj(eigValArr_D2);

applyD1 = @(x) applyCyclicConv2D(x,eigValArr_D1);
applyD1Trans = @(x) applyCyclicConv2D(x,eigValArr_D1Trans);
applyD2 = @(x) applyCyclicConv2D(x,eigValArr_D2);
applyD2Trans = @(x) applyCyclicConv2D(x,eigValArr_D2Trans);

applyD = @(x) cat(3,applyD1(x),applyD2(x));
applyDTrans = @(y) applyD1Trans(y(:,:,1)) + applyD2Trans(y(:,:,2));

[R1,R2] = makeSparseGradRemainderMatrices(numRows,numCols);
applyR1 = @(x) reshape(R1*x(:),numRows,numCols);
applyR1Trans = @(x) reshape(R1'*x(:),numRows,numCols);
applyR2 = @(x) reshape(R2*x(:),numRows,numCols);
applyR2Trans = @(x) reshape(R2'*x(:),numRows,numCols);

applyR = @(x) cat(3,applyR1(x),applyR2(x));
applyRTrans = @(y) applyR1Trans(y(:,:,1)) + applyR2Trans(y(:,:,2));

applyG = @(x) applyD(x) + applyR(x);
applyGTrans = @(y) applyDTrans(y) + applyRTrans(y);

applyK = @(x) cat(3,applyA(x),applyG(x));
applyKTrans = @(z) applyATrans(z(:,:,1)) + applyGTrans(z(:,:,2:3));

%%%%%%%%%%%%%%%%%%%

nrmA = max(abs(eigValArr_P(:)));
nrmK = sqrt(8 + nrmA^2);

% nrmK_est = estimateNormByPowerIteration(applyK,applyKTrans,b,500);

tau = params.tau;
sigma = 1/(tau*nrmK^2);

xn = b;
yn = applyK(xn);
xnBar = xn;

costs = [];
figure('Name','xn inside Chambolle-Pock')

for n = 1:maxIter 
    
    ATransyn_1 = applyATrans(yn(:,:,1));
    GTransyn_2 = applyGTrans(yn(:,:,2:end));
    KTransyn = ATransyn_1 + GTransyn_2;          
    
    % compute proxTauG, which simply projects onto C.
    xnp1 = min(xn - tau*KTransyn,1);
    xnp1 = max(xnp1,0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
            
    z = yn + sigma*applyK(2*xnp1 - xn);
    
    % compute proxSigmaFStar evaluated at z.
    z1 = z(:,:,1);
    z2 = z(:,:,2:3);
    
    term = z1 - sigma*b;    
    ynp1_1 = min(term,1);
    ynp1_1 = max(ynp1_1,-1);
        
    ynp1_2 = gamma*projectOntoDualIsoNormUnitBall(z2/gamma);
    ynp1 = cat(3,ynp1_1,ynp1_2); % This is dual feasible.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    xn = overRelax*xnp1 + (1-overRelax)*xn;
    yn = overRelax*ynp1 + (1-overRelax)*yn; 
        
    Axnp1 = applyA(xnp1);
    Gxnp1 = applyG(xnp1);    
            
    cost = sum(sum(abs(Axnp1-b))) + gamma*evalIsoNorm(Gxnp1);
    costs = [costs,cost];    
    
    if mod(n,showTrigger) == 0
        imshow(xn,[])
        disp(['Chambolle-Pock iteration is: ',num2str(n)])
        disp(['primal cost is: ',num2str(cost)])         
        keyboard        
    end     
    
end

end