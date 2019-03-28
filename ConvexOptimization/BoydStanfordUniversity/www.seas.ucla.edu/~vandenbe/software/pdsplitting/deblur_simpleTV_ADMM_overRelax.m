function [xk,costs] = deblur_simpleTV_ADMM_overRelax(b,mask,params)

% We solve minimize rho*||Ax - b||_1 + || Dx ||_{iso}
% subject to 0 <= x <= 1.
% A convolves with mask using periodic boundary conditions.
% D is a discrete gradient operator.

rho = params.rho;
t = params.stepSize;
maxIter = params.maxIter;
showTrigger = params.showTrigger;
beta = params.beta; 
alpha = params.alpha; % over-relaxation parameter

%%% Using two stepsizes is equivalent to modifying objective and the
%%% operators, and using a single step size t.
%%% K and D are multiplied by beta, and mu and b 
%%% are adjusted to compensate for that.
mask = mask*beta;
rho = rho/beta;
b = b*beta;
%%% Finished modifying objective.
% After this modification, we are now minimizing
% rho*||Ax - b||_1 + (1/beta)*||Dx||_{iso} s.t. 0 <= x <= 1  (using the notation of this code)
% We do this using Douglas-Rachford with a single step size t.

evalNorm = @(y) sum(sum(sqrt(y(:,:,1).^2 + y(:,:,2).^2)));

[numRows,numCols] = size(b);

adjMask = fliplr(mask);
adjMask = flipud(adjMask);

eigValArr_A = eigValArrForCyclicConvOp(mask,numRows,numCols);
eigValArr_D1 = eigValArrForCyclicConvOp(beta*[-1 1]',numRows,numCols); % Note that D is being multiplied by beta
eigValArr_D2 = eigValArrForCyclicConvOp(beta*[-1 1],numRows,numCols); % as mentioned above.

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

applyC = @(x) cat(3,applyA(x),applyD(x));
applyCTrans = @(z) applyATrans(z(:,:,1)) + applyDTrans(z(:,:,2:3));

applyMtrx = @(x) applyATrans(applyA(x)) + applyDTrans(applyD(x)) + x; 
eigValsMtrx = eigValArr_ATrans.*eigValArr_A ...
                + eigValArr_D1Trans.*eigValArr_D1 + eigValArr_D2Trans.*eigValArr_D2 + ones(numRows,numCols);
            
xk = b;
yk = applyD(xk);
uk = xk;
vk = applyA(xk) - b;

zHat = zeros(size(yk));
pHat = zeros(size(xk));
qHat = zeros(size(b));

costs = [];

figure('Name','xk inside ADMM')

for k = 1:maxIter
    
    rhs = uk - pHat/t + applyATrans(b + vk - qHat/t) + applyDTrans(yk - zHat/t);
    xkp1 = ifft2(fft2(rhs)./eigValsMtrx);
    
%     check = applyMtrx(xkp1) - rhs;
%     mx = max(abs(check(:)));

    Axkp1 = applyA(xkp1);
    Dxkp1 = applyD(xkp1);
       
    ykp1 = proxIsoNorm(1/(beta*t),alpha*Dxkp1 + (1-alpha)*yk + zHat/t);
%     ykp1 = proxIsoNorm(1/(beta*t),applyD(xkp1) + zHat/t);     
        
    ukp1 = min(alpha*xkp1 + (1-alpha)*uk + pHat/t,1);    
%     ukp1 = min(xkp1 + pHat/t,1);    
    ukp1 = max(ukp1,0);
    
    vkp1 = prox1Norm(alpha*Axkp1 + (1-alpha)*(vk + b) - b + qHat/t,rho/t);
%     vkp1 = prox1Norm(applyA(xkp1) - b + qHat/t,rho/t);

    pHat = pHat + t*(alpha*xkp1 + (1-alpha)*uk - ukp1);
    qHat = qHat + t*(alpha*Axkp1 + (1-alpha)*(vk + b) - b - vkp1);
    zHat = zHat + t*(alpha*Dxkp1 + (1-alpha)*yk - ykp1);
    
    xk = xkp1;
    yk = ykp1;
    uk = ukp1;
    vk = vkp1;
    
%     zHat = zHat + t*(applyD(xk) - yk);
%     pHat = pHat + t*(xk - uk);
%     qHat = qHat + t*(applyA(xk) - b - vk);
            
    AukMinusb = applyA(uk) - b;
    Duk = applyD(uk);
    primalCost = rho*sum(abs(AukMinusb(:))) + (1/beta)*evalNorm(Duk);           
    costs = [costs,primalCost];      
    
    if mod(k,showTrigger) == 0 
        imshow(xk,[])
        disp(['ADMM iteration is: ',num2str(k)])
        disp(['primal cost is: ',num2str(primalCost)])         
        keyboard        
    end                  

end