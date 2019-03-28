function [xPlus_x1,costs] = deblur_mixedTightFrame_drPrimal(b,kernel,Ks,params)
               
% This code corresponds to section 6.2 in the paper
% "Image deblurring by primal-dual operator splitting." 
% We solve minimize .5*||Kx - b||_F^2 + gamma*|| Dx ||_1.
% K = Kp + Ks, Ks is sparse, Kp uses periodic BCs.
% D is a tight frame analysis operator.

% In this code primal and dual step sizes are handled as follows:
% the operator A gets multiplied by beta, and g is replaced by
% gTilde, as in p. 7 of the paper.  
% With these modifications to the objective function,
% we then use Douglas-Rachford with the single step size t.

% After modifying K,D, and b, our problem is to 
% minimize (.5/beta^2)||Kx - b||_^2 + (gamma/beta)||Dx||_1.

maxIter = params.maxIter;
t = params.t; % primal step size
beta = params.beta;
gamma = params.gamma;
showTrigger = params.showTrigger;
overRelax = params.overRelax; 

[numRows,numCols] = size(b);
numPix = numRows*numCols;

Ks = beta*Ks;
kernel = beta*kernel;
b = beta*b;

eigValArr_Kp = eigValArrForCyclicConvOp(kernel,numRows,numCols);
eigValArr_KpTrans = conj(eigValArr_Kp);

applyKp = @(x) applyCyclicConv2D(x,eigValArr_Kp);
applyKpTrans = @(x) applyCyclicConv2D(x,eigValArr_KpTrans);

applyKs = @(x) reshape(Ks*x(:),numRows,numCols);
applyKsTrans = @(x) reshape(Ks'*x(:),numRows,numCols);

applyK = @(x) applyKp(x) + applyKs(x);
applyKTrans = @(x) applyKpTrans(x) + applyKsTrans(x);

%%%%%
scale = [3 3 3 4 4]; 
qmf1 = MakeONFilter('Symmlet',4);
qmf2 = MakeONFilter('Symmlet',4);
ndir = 0;
alpha = 2^(ndir+2)+2; % D^T D = alpha I.
alpha = (beta^2)*alpha; % This is because we have TWO STEPSIZES.
                        % When we modify D, alpha is modified also.                                                

applyD = @(x) beta*shearletTransform(x,qmf1,qmf2,scale,ndir);
applyDTrans = @(y) beta*shearletAdjTransform(y,qmf1,qmf2,scale,ndir);

applyB = @(x) cat(3,applyKp(x),applyD(x));
applyBTrans = @(y) applyKpTrans(y(:,:,1)) + applyDTrans(y(:,:,2:end));

dim = size(applyD(b),3);
applyC = @(x) cat(3,applyKs(x),zeros(numRows,numCols,dim));
applyCTrans = @(y) applyKsTrans(y(:,:,1));

applyA = @(x) cat(3,applyK(x),applyD(x));
applyATrans = @(y) applyKTrans(y(:,:,1)) + applyDTrans(y(:,:,2:end));

applyMtrx = @(x) (1 + alpha)*x + applyKpTrans(applyKp(x));
eigValsMtrx = eigValArr_KpTrans.*eigValArr_Kp + (1+alpha)*ones(numRows,numCols);
            
sparseMtrx = speye(numPix) + Ks'*Ks ;

[U,statusInfo,perm] = chol(sparseMtrx,'vector');
UTrans = U';
permInv(perm) = 1:numPix;

p_x1 = b;
p_x2 = b;
p_y1 = applyB(b);
p_y2 = applyC(b);

costs = [];
figure('Name','xPlus_x inside DR')

for n = 1:maxIter 

    % Evaluate resolvent at p_x1,p_x2,p_y1,p_y2.    
    xPlus_x1 = (p_x1 + p_x2)/2;
    xPlus_x2 = xPlus_x1;
    
    term = p_y1 + p_y2; % z is prox op. of g evaluated at term, with parameter 2t.
    termData = term(:,:,1);
    termReg = term(:,:,2:end);
    
%     zData = prox1Norm(termData - b,2*t/beta) + b;
    zData = (b/beta^2 + termData/(2*t))/(1/beta^2 + 1/(2*t));
    zReg = prox1Norm(termReg,2*t*gamma/beta);        
    z = cat(3,zData,zReg);    
    xPlus_y1 = (z + p_y1 - p_y2)/2;
    xPlus_y2 = (z + p_y2 - p_y1)/2;
    
    x1Hat = 2*xPlus_x1 - p_x1;
    x2Hat = 2*xPlus_x2 - p_x2;
    y1Hat = 2*xPlus_y1 - p_y1;
    y2Hat = 2*xPlus_y2 - p_y2;
    
    rhs = x1Hat + applyBTrans(y1Hat);
    yPlus_x1 = ifft2(fft2(rhs)./eigValsMtrx);
    yPlus_y1 = applyB(yPlus_x1);
    
    rhs2 = x2Hat + applyCTrans(y2Hat);
    rhs2 = rhs2(:); 
    temp = UTrans\rhs2(perm);    
    w = U\temp;    
    yPlus_x2 = w(permInv);
    yPlus_x2 = reshape(yPlus_x2,[numRows,numCols]);   
    yPlus_y2 = applyC(yPlus_x2);
    
    check = sparseMtrx*yPlus_x2(:) - rhs2(:);
    check = max(abs(check(:)));
    
    pPlus_x1 = p_x1 + overRelax*(yPlus_x1 - xPlus_x1);
    pPlus_x2 = p_x2 + overRelax*(yPlus_x2 - xPlus_x2);
    pPlus_y1 = p_y1 + overRelax*(yPlus_y1 - xPlus_y1);
    pPlus_y2 = p_y2 + overRelax*(yPlus_y2 - xPlus_y2);
    
    p_x1 = pPlus_x1;
    p_x2 = pPlus_x2;
    p_y1 = pPlus_y1;
    p_y2 = pPlus_y2;
               
    Dx = applyD(xPlus_x1);
    KxMinusb = applyK(xPlus_x1) - b;    
    primalCost = (.5/beta^2)*sum(KxMinusb(:).^2) + (gamma/beta)*sum(abs(Dx(:)));
    primalCost = primalCost/gamma; % this is to compare with CP costs.
    
    costs = [costs,primalCost];         
                    
    if mod(n,showTrigger) == 0         
        imshow(xPlus_x1,[])    
        disp(['drPrimal iteration is: ',num2str(n)])
        disp(['primal cost is: ',num2str(primalCost)])         
        keyboard        
    end                 
                
end


end