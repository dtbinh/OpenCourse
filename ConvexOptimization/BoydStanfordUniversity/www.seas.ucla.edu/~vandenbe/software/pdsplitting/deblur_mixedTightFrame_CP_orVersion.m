function [xn,costs] = deblur_mixedTightFrame_CP_orVersion(b,mask,E,params)
% orVersion stands for "overrelaxed version" of Chambolle-Pock.

% We solve minimize .5*rho*||Ax - b||^2 + || Dx ||_1
% A = B + E, B convolves with mask using periodic BCs, E is sparse.
% D is analysis operator for shearlet tight frame.

maxIter = params.maxIter;
showTrigger = params.showTrigger;
rho = params.rho;
overRelax = params.overRelax;

[numRows,numCols] = size(b);
%%%%%%%%%%%% Set up operators.
adjMask = fliplr(mask);
adjMask = flipud(adjMask);

eigValArr_B = eigValArrForCyclicConvOp(mask,numRows,numCols);
eigValArr_BTrans = conj(eigValArr_B);

applyB = @(x) applyCyclicConv2D(x,eigValArr_B);
applyBTrans = @(x) applyCyclicConv2D(x,eigValArr_BTrans);

applyE = @(x) reshape(E*x(:),numRows,numCols);
applyETrans = @(x) reshape(E'*x(:),numRows,numCols);

applyA = @(x) applyB(x) + applyE(x);
applyATrans = @(x) applyBTrans(x) + applyETrans(x);

scale = [3 3 3 4 4]; 
qmf1 = MakeONFilter('Symmlet',4);
qmf2 = MakeONFilter('Symmlet',4);
ndir = 0;
alpha = 2^(ndir+2)+2; % D^T D = alpha I.

applyD = @(x) shearletTransform(x,qmf1,qmf2,scale,ndir);
applyDTrans = @(y) shearletAdjTransform(y,qmf1,qmf2,scale,ndir);

applyK = @(x) cat(3,applyA(x),applyD(x));
applyKTrans = @(z) applyATrans(z(:,:,1)) + applyDTrans(z(:,:,2:end));

%%%%%%%%%%%%%%%%%%%

nrmB = max(abs(eigValArr_B(:)));
nrmK = sqrt(alpha + nrmB^2);

% [nrmK_est,powerIters] = estimateNormByPowerIteration(applyK,applyKTrans,b,500);

tau = params.tau;
sigma = 1/(tau*nrmK^2);

xn = b;
yn = applyK(xn);
xnBar = xn;

costs = [];
figure('Name','xn inside Chambolle-Pock')

for n = 1:maxIter

	KTransyn = applyKTrans(yn);       
    xnp1 = xn - tau*KTransyn;
            
    z = yn + sigma*applyK(2*xnp1 - xn);
    
    % compute proxSigmaFStar evaluated at z.
    z1 = z(:,:,1);
    z2 = z(:,:,2:end);
    
    ynp1_1 = (rho*z1 - (rho*sigma)*b)/(rho + sigma);    
        
    ynp1_2 = min(z2,1);
    ynp1_2 = max(ynp1_2,-1);
    ynp1 = cat(3,ynp1_1,ynp1_2);   
    
	xn = overRelax*xnp1 + (1-overRelax)*xn;
	yn = overRelax*ynp1 + (1-overRelax)*yn;
        
    Dxnp1 = applyD(xnp1);
    Axnp1Minusb = applyA(xnp1) - b;
    cost = .5*rho*sum(Axnp1Minusb(:).^2) + sum(abs(Dxnp1(:))); 
    costs = [costs,cost];                               
    
    if mod(n,showTrigger) == 0
        imshow(xn,[])
        disp(['Chambolle-Pock iteration is: ',num2str(n)])
        disp(['primal cost is: ',num2str(cost)])         
        keyboard         
    end      
    
end

end