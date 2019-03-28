function [xn,costs] = deblur_simpleTV_CP_orVersion(b,mask,params)

% We solve minimize rho*||Ax - b||_1 + || Dx ||_{iso}
% subject to 0 <= x <= 1.
% A convolves with mask using periodic boundary conditions.
% D is a discrete gradient operator.

% F(y1,y2) = rho*||y1 - b||_1 + || y2 ||_{iso}.
% F1(y1) = rho*||y1 - b||_1 and F2(y2) = || y2 ||_{iso}.
% F1Star(z1) = deltaB( z1/rho) + < b , z1 >
% where B is the unit ball for the infinity norm.
% deltaB is indicator function for B.
% prox_{sigma F1Star(z1Hat)} = projection of z1Hat - sigma*b onto rho*B.

maxIter = params.maxIter;
showTrigger = params.showTrigger;
rho = params.rho;
overRelax = params.overRelax;

evalNorm = @(y) sum(sum(sqrt(y(:,:,1).^2 + y(:,:,2).^2)));

[numRows,numCols] = size(b);
%%%%%%%%%%%% Set up operators.
eigValArr_A = eigValArrForCyclicConvOp(mask,numRows,numCols);
eigValArr_D1 = eigValArrForCyclicConvOp([-1 1]',numRows,numCols);
eigValArr_D2 = eigValArrForCyclicConvOp([-1 1],numRows,numCols);

eigValArr_ATrans = conj(eigValArr_A);
eigValArr_D1Trans = conj(eigValArr_D1);
eigValArr_D2Trans = conj(eigValArr_D2);

applyA = @(x) applyCyclicConv2D(x,eigValArr_A);
applyATrans = @(x) applyCyclicConv2D(x,eigValArr_ATrans);
applyD1 = @(x) applyCyclicConv2D(x,eigValArr_D1);
applyD1Trans = @(x) applyCyclicConv2D(x,eigValArr_D1Trans);
applyD2 = @(x) applyCyclicConv2D(x,eigValArr_D2);
applyD2Trans = @(x) applyCyclicConv2D(x,eigValArr_D2Trans);

applyD = @(x) cat(3,applyD1(x),applyD2(x));
applyDTrans = @(y) applyD1Trans(y(:,:,1)) + applyD2Trans(y(:,:,2));

applyK = @(x) cat(3,applyA(x),applyD(x));
applyKTrans = @(z) applyATrans(z(:,:,1)) + applyDTrans(z(:,:,2:3));

%%%%%%%%%%%%%%%%%%%

nrmA = max(abs(eigValArr_A(:)));
nrmK = sqrt(8 + nrmA^2);

% nrmK = estimateNormByPowerIteration(applyK,applyKTrans,b,500);
% disp(['nrmK (estimated by power iteration) is: ',num2str(nrmK)]) % 2.8269

tau = params.tau;
sigma = 1/(tau*nrmK^2);

xn = b;
yn = applyK(xn);

costs = [];
figure('Name','xn inside Chambolle-Pock')

for n = 1:maxIter 

    ATransyn_1 = applyATrans(yn(:,:,1));
    DTransyn_2 = applyDTrans(yn(:,:,2:end));
    KTransyn = ATransyn_1 + DTransyn_2;
                            
    xnp1 = min(xn - tau*KTransyn,1);
    xnp1 = max(xnp1,0);                                     
            
    z = yn + sigma*applyK(2*xnp1 - xn);
    
    % compute proxSigmaFStar evaluated at z.
    z1 = z(:,:,1);
    z2 = z(:,:,2:3);        
            
    term = z1 - sigma*b;
    ynp1_1 = min(term,rho);
    ynp1_1 = max(ynp1_1,-rho); 
    
    ynp1_2 = projectOntoDualIsoNormUnitBall(z2);    
        
    ynp1 = cat(3,ynp1_1,ynp1_2);   
    
	xn = overRelax*xnp1 + (1-overRelax)*xn;
	yn = overRelax*ynp1 + (1-overRelax)*yn;	
         
    
    Axnp1 = applyA(xnp1);
    Dxnp1 = applyD(xnp1);          
        
    cost = rho*sum(sum(abs(Axnp1 - b))) + evalNorm(Dxnp1);
    costs = [costs,cost];                     
                 
    if mod(n,showTrigger) == 0       
        imshow(xn,[])
        disp(['Chambolle-Pock iteration is: ',num2str(n)])
        disp(['primal cost is: ',num2str(cost)])         
        keyboard
    end        
    
end

end