function [xk,costs] = deblur_simpleTightFrame_ADMM(b,mask,params)
% We solve minimize rho*||Ax - b||_1 + || Dx ||_1
% subject to 0 <= x <= 1.
% A convolves with mask using periodic boundary conditions.
% D is the analysis operator for a shearlet tight frame.

t = params.stepSize;
maxIter = params.maxIter;
showTrigger = params.showTrigger;
rho = params.rho;

[numRows,numCols] = size(b);

adjMask = fliplr(mask);
adjMask = flipud(adjMask);

eigValArr_A = eigValArrForCyclicConvOp(mask,numRows,numCols);
eigValArr_ATrans = conj(eigValArr_A);

applyA = @(x) applyCyclicConv2D(x,eigValArr_A);
applyATrans = @(x) applyCyclicConv2D(x,eigValArr_ATrans);

scale = [3 3 3 4 4]; 
qmf1 = MakeONFilter('Symmlet',4);
qmf2 = MakeONFilter('Symmlet',4);
ndir = 0;
alpha = 2^(ndir+2)+2; % D^T D = alpha I.

applyD = @(x) shearletTransform(x,qmf1,qmf2,scale,ndir);
applyDTrans = @(y) shearletAdjTransform(y,qmf1,qmf2,scale,ndir);

applyC = @(x) cat(3,applyA(x),applyD(x));
applyCTrans = @(z) applyATrans(z(:,:,1)) + applyDTrans(z(:,:,2:end));

applyMtrx = @(x) applyATrans(applyA(x)) + (alpha+1)*x; 
eigValsMtrx = eigValArr_ATrans.*eigValArr_A + (alpha+1)*ones(numRows,numCols);
            
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
        
    ykp1 = prox1Norm(applyD(xkp1) + zHat/t,1/t);
    
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
   
    Duk = applyD(uk);
    AukMinusb = applyA(uk) - b;
    primalCost = rho*sum(abs(AukMinusb(:))) + sum(abs(Duk(:)));        
    costs = [costs,primalCost];      
        
    if mod(k,showTrigger) == 0        
        imshow(xk,[]) 
        disp(['ADMM iteration is: ',num2str(k)])
        disp(['primal cost is: ',num2str(primalCost)])         
        keyboard        
    end    
        
end