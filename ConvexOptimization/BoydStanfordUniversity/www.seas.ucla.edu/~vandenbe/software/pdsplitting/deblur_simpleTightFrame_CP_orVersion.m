function [xn,costs] = deblur_simpleTightFrame_CP_orVersion(b,mask,params)
% orVersion stands for "overrelaxed version" of Chambolle-Pock.

% We solve minimize rho*||Ax - b||_1 + || Dx ||_1
% subject to 0 <= x <= 1.
% A convolves with mask using periodic boundary conditions.
% D is analysis operator for a shearlet tight frame.

% F(y1,y2) = rho*||y1 - b||_1 + || y2 ||_1.
% F1(y1) = rho*||y1 - b||_1 and F2(y2) = || y2 ||_1.
% F1Star(z1) = deltaB( z1/rho) + < b , z1 >
% where B is the unit ball for the infinity norm.
% deltaB is indicator function for B.
% prox_{sigma F1Star(z1Hat)} = projection of z1Hat - sigma*b onto rho*B.

% F2Star(z2) = deltaB(z2);
% prox_{sigma F2Star(z2Hat)} = projection of z2Hat onto B.

maxIter = params.maxIter;
showTrigger = params.showTrigger;
rho = params.rho;
overRelax = params.overRelax;

[numRows,numCols] = size(b);

%%%%%%%%%%%% Set up operators.
eigValArr_A = eigValArrForCyclicConvOp(mask,numRows,numCols);
eigValArr_ATrans = conj(eigValArr_A);

applyA = @(x) applyCyclicConv2D(x,eigValArr_A);
applyATrans = @(x) applyCyclicConv2D(x,eigValArr_ATrans);

% set up D using Shearlab.
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

nrmA = max(abs(eigValArr_A(:)));
nrmK = sqrt(alpha + nrmA^2);

% [nrmK_est,powerIters] = estimateNormByPowerIteration(applyK,applyKTrans,b,500);

tau = params.tau;
sigma = 1/(tau*nrmK^2);

xn = b;
yn = applyK(xn);

costs = [];
figure('Name','xn inside CP')

for n = 1:maxIter 

    KTransyn = applyKTrans(yn);
                 
    xnp1 = min(xn - tau*KTransyn,1);
    xnp1 = max(xnp1,0);      %%%%%                               
            
    z = yn + sigma*applyK(2*xnp1 - xn);
    
    % compute proxSigmaFStar evaluated at z.
    z1 = z(:,:,1);
    z2 = z(:,:,2:end);            
        
    term = z1 - sigma*b;
    ynp1_1 = min(term,rho);
    ynp1_1 = max(ynp1_1,-rho);     
    
    ynp1_2 = min(z2,1);
    ynp1_2 = max(ynp1_2,-1);
    ynp1 = cat(3,ynp1_1,ynp1_2);  
    
	xn = overRelax*xnp1 + (1-overRelax)*xn;
	yn = overRelax*ynp1 + (1-overRelax)*yn;	    
    
    Axnp1 = applyA(xnp1);
    Dxnp1 = applyD(xnp1);    
       
    cost = rho*sum(sum(abs(Axnp1 - b))) + sum(abs(Dxnp1(:)));    
    costs = [costs,cost];                     
                       
    if mod(n,showTrigger) == 0
        imshow(xn,[])    
        disp(['Chambolle-Pock iteration is: ',num2str(n)])
        disp(['primal cost is: ',num2str(cost)])         
        keyboard        
    end     
    
end

end