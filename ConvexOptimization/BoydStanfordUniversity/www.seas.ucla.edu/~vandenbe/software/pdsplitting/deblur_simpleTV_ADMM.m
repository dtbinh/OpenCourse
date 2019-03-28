function [xk,costs] = deblur_simpleTV_ADMM(b,mask,params)

% We solve minimize rho*||Ax - b||_1 + || Dx ||_{iso}
% subject to 0 <= x <= 1.
% A convolves with mask using periodic boundary conditions.
% D is a discrete gradient operator.

rho = params.rho;
t = params.stepSize;
maxIter = params.maxIter;
showTrigger = params.showTrigger;

evalNorm = @(y) sum(sum(sqrt(y(:,:,1).^2 + y(:,:,2).^2)));

[numRows,numCols] = size(b);

adjMask = fliplr(mask);
adjMask = flipud(adjMask);

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
    
    ykp1 = proxIsoNorm(1/t,applyD(xkp1) + zHat/t);    
        
    ukp1 = min(xkp1 + pHat/t,1);
    ukp1 = max(ukp1,0);
    
    vkp1 = prox1Norm(applyA(xkp1) - b + qHat/t,rho/t);
    
    xk = xkp1;
    yk = ykp1;
    uk = ukp1;
    vk = vkp1;
    
    zHat = zHat + t*(applyD(xk) - yk);
    pHat = pHat + t*(xk - uk);
    qHat = qHat + t*(applyA(xk) - b - vk);
            
    AukMinusb = applyA(uk) - b;
    Duk = applyD(uk);
    primalCost = rho*sum(abs(AukMinusb(:))) + evalNorm(Duk);           
    costs = [costs,primalCost];      
    
    if mod(k,showTrigger) == 0 
        imshow(xk,[])
        disp(['ADMM iteration is: ',num2str(k)])
        disp(['primal cost is: ',num2str(primalCost)])         
        keyboard        
    end                  

end