function [xk,costs] = deblur_simpleTightFrame_ADMM_overRelax(b,mask,params)

% We solve minimize rho*||Ax - b||_1 + || Dx ||_1
% subject to 0 <= x <= 1.
% A convolves with mask using periodic boundary conditions.
% D is a discrete gradient operator.

rho = params.rho;
t = params.stepSize;
maxIter = params.maxIter;
showTrigger = params.showTrigger;
beta = params.beta; 
gamma = params.gamma; % over-relaxation parameter

%%% Using two stepsizes is equivalent to modifying objective and the
%%% operators, and using a single step size t.
%%% K and D are multiplied by beta, and mu and b 
%%% are adjusted to compensate for that.
mask = mask*beta;
rho = rho/beta;
b = b*beta;
%%% Finished modifying objective.
% After this modification, we are now minimizing
% rho*||Ax - b||_1 + (1/beta)*||Dx||_1 s.t. 0 <= x <= 1  (using the notation of this code)
% We do this using Douglas-Rachford with a single step size t.

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
alpha = (beta^2)*alpha; % This is because we have TWO STEPSIZES.
                        % When we modify D, alpha is modified also.

applyD = @(x) beta*shearletTransform(x,qmf1,qmf2,scale,ndir);
applyDTrans = @(y) beta*shearletAdjTransform(y,qmf1,qmf2,scale,ndir);

applyMtrx = @(x) applyATrans(applyA(x)) + (alpha+1)*x; 
eigValsMtrx = eigValArr_ATrans.*eigValArr_A ...
                + (alpha+1)*ones(numRows,numCols);
            
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
           
    ykp1 = prox1Norm(gamma*Dxkp1 + (1-gamma)*yk + zHat/t,1/(beta*t));  
        
    ukp1 = min(gamma*xkp1 + (1-gamma)*uk + pHat/t,1);    
    ukp1 = max(ukp1,0);
    
    vkp1 = prox1Norm(gamma*Axkp1 + (1-gamma)*(vk + b) - b + qHat/t,rho/t);

    pHat = pHat + t*(gamma*xkp1 + (1-gamma)*uk - ukp1);
    qHat = qHat + t*(gamma*Axkp1 + (1-gamma)*(vk + b) - b - vkp1);
    zHat = zHat + t*(gamma*Dxkp1 + (1-gamma)*yk - ykp1);
    
    xk = xkp1;
    yk = ykp1;
    uk = ukp1;
    vk = vkp1;    
            
    AukMinusb = applyA(uk) - b;
    Duk = applyD(uk);
    primalCost = rho*sum(abs(AukMinusb(:))) + (1/beta)*sum(abs(Duk(:)));           
    costs = [costs,primalCost];      
    
    if mod(k,showTrigger) == 0 
        imshow(xk,[])
        disp(['ADMM iteration is: ',num2str(k)])
        disp(['primal cost is: ',num2str(primalCost)])         
        keyboard        
    end                  

end