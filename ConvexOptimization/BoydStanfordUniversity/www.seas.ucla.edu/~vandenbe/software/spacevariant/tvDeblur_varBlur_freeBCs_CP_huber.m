function [xn,costs,errors] = tvDeblur_varBlur_freeBCs_CP_huber(b,kernels,U,params,xStar)

% We minimize phif(U1*K1 + ... + U_P*K_P)x - b) + gamma*||Dx||_{iso}.
% phif = \sum_i phifi.
% phifi(u) = huber(u - bi) unless i is a border pixel.  In that case phifi = 0.
% Objective is F(Kx) + G(x) where G(x) = 0.
% A = U_1*K_1 + ... + U_P*K_P.
% Each Kp convolves using periodic BCs.
% K = [ A ; D].  (Not the best notation because K doesn't relate to
% K1,...,K_P as expected.)

maxIter = params.maxIter;
showTrigger = params.showTrigger;
gamma = params.gamma;
overRelax = params.overRelax;
mu = params.mu;  % Huber penalty parameter.

[numRows,numCols,numKernels] = size(U);

d = (size(kernels,1) - 1)/2 ;
evalIsoNorm = @(y) sum(sum(sqrt(y(:,:,1).^2 + y(:,:,2).^2)));

%%%%%%%
Ureduced = U;
U = (1/numKernels)*ones(numRows + 2*d,numCols + 2*d,numKernels);
for k = 1:numKernels
    U(d+1:d+numRows,d+1:d+numCols,k) = Ureduced(:,:,k);
end
%%%%%%%

eigValArrs_Kp = zeros(numRows + 2*d,numCols + 2*d,numKernels);
eigValArrs_KpTrans = zeros(numRows + 2*d,numCols + 2*d,numKernels);
for k = 1:numKernels
    eigValArr = eigValArrForCyclicConvOp(kernels(:,:,k),numRows + 2*d,numCols + 2*d);
    eigValArrs_Kp(:,:,k) = eigValArr;
    eigValArrs_KpTrans(:,:,k) = conj(eigValArr);
end

%%%%%%%%%%%% Set up operators.
% Create applyA
applyA = @(x) applyNagyBlur(x,eigValArrs_Kp,U)
applyATrans = @(y) applyNagyBlurTrans(y,eigValArrs_KpTrans,U);

% Create applyD
eigValArr_D1 = eigValArrForCyclicConvOp([-1 1]',numRows + 2*d,numCols + 2*d); 
eigValArr_D2 = eigValArrForCyclicConvOp([-1 1],numRows + 2*d,numCols + 2*d); 

eigValArr_D1Trans = conj(eigValArr_D1);
eigValArr_D2Trans = conj(eigValArr_D2);

applyD1 = @(x) applyCyclicConv2D(x,eigValArr_D1);
applyD1Trans = @(x) applyCyclicConv2D(x,eigValArr_D1Trans);
applyD2 = @(x) applyCyclicConv2D(x,eigValArr_D2);
applyD2Trans = @(x) applyCyclicConv2D(x,eigValArr_D2Trans);

applyD = @(x) cat(3,applyD1(x),applyD2(x));
applyDTrans = @(y) applyD1Trans(y(:,:,1)) + applyD2Trans(y(:,:,2));

% Create applyK
applyK = @(x) cat(3,applyA(x),applyD(x));
applyKTrans = @(z) applyATrans(z(:,:,1)) + applyDTrans(z(:,:,2:end));

%%%%%%%%%%%%%%%%%%%

bHat = zeros(numRows + 2*d,numCols + 2*d);
bHat(d+1:d+numRows,d+1:d+numCols) = b;

% either compute or load a precomputed nrmK
if 0
    [nrmK_est,powerIters] = estimateNormByPowerIteration(applyK,applyKTrans,bHat,2000);
else
    load 'nrmK_CP_simple.mat'
    disp('LOADING a precomputed nrmK !')    
end

nrmK = nrmK_est;

tau = params.tau;
sigma = 1/(tau*nrmK^2);

xn = bHat;
yn = applyK(xn);

costs = [];
errors = [];

for n = 1:maxIter

	KTransyn = applyKTrans(yn);       
    xnp1 = xn - tau*KTransyn; % prox operator of G is just identity operator, since G = 0.
            
    z = yn + sigma*applyK(2*xnp1 - xn);
    
    % compute proxSigmaFStar evaluated at z. 
    z1 = z(:,:,1);
    z2 = z(:,:,2:end);
    
    %%%%%%%  Compute proxg1Star at z1.    
    proxg1 = z1/sigma;
    in = proxg1(d+1:d+numRows,d+1:d+numCols);
    proxg1(d+1:d+numRows,d+1:d+numCols) = b + proxHuber(in - b,mu,1/sigma);
    proxg1Star = z1 - sigma*proxg1;
    
    % Next three lines were for the version with L2 data fidelity term.    
%     proxg1 = z1/sigma;
%     proxg1(d+1:d+numRows,d+1:d+numCols) = (b + z1(d+1:d+numRows,d+1:d+numCols))/(sigma + 1);    
%     
%     proxg1Star = z1 - sigma*proxg1;        

    %%%%%%% Finished computing proxg1Star at z1.    
    
    proxg2Star = gamma*projectOntoDualIsoNormUnitBall(z2/gamma); 
    
    ynp1_1 = proxg1Star;
    ynp1_2 = proxg2Star;
    ynp1 = cat(3,ynp1_1,ynp1_2);
        
    %%%%%%         
	xn = overRelax*xnp1 + (1-overRelax)*xn;
	yn = overRelax*ynp1 + (1-overRelax)*yn;           

    % Compute objective function value and error.    
       
    Ax = applyA(xn);    
    Dx = applyD(xn);  
    
    in = Ax(d+1:d+numRows,d+1:d+numCols) - b;    
    prx = prox1Norm(in, mu);
    huberTerm = abs(prx) + (.5/mu)*(prx - in).^2;
    huberTerm = sum(huberTerm(:));
    cost = huberTerm + gamma*evalIsoNorm(Dx); % obj. fun. value
    costs = [costs,cost];
    
    if nargin == 5
        error = xn - xStar; 
        error = norm(error(:));            
        errors = [errors,error];
    end
    

    if mod(n,showTrigger) == 0
        
        if n == showTrigger
            figure('Name','xn inside Chambolle-Pock')
        end
        
        reduced = xn(d+1:d+numRows,d+1:d+numCols);
        
        imshow(reduced,[])
        disp(['Chambolle-Pock iteration is: ',num2str(n)])
        disp(['primal cost is: ',num2str(cost)])         
%         keyboard         
    end      
    
end

end