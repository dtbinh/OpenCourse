function [uk,costs] = deblur_mixedTightFrame_admm(b,kernel,Ks,params)
               
% This code corresponds to section 5.2 in the paper
% "Image deblurring by primal-dual operator splitting." 
% We solve minimize .5*|| Kx - b||_2^2 + gamma*|| Dx ||_1.
% K = Kp + Ks. 
% D is the analysis operator for a shearlet tight frame.

% In this code primal and dual step sizes are handled as follows:
% the operator A gets multiplied by beta, and g is replaced by
% gTilde, as in p. 7 of the paper.  
% With these modifications to the objective function,
% we then use Douglas-Rachford with the single step size t.

% After modifying K,D, and b, our problem is to 
% minimize (.5/beta^2)||Kx - b||_2^2 + (gamma/beta)||Dx||_1.

maxIter = params.maxIter;
t = params.t; % primal step size
beta = params.beta; % for scaling the problem.
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
eigValsMtrx = eigValArr_KpTrans.*eigValArr_Kp + (1 + alpha)*ones(numRows,numCols);
            
sparseMtrx = speye(numPix) + Ks'*Ks;

[U,statusInfo,perm] = chol(sparseMtrx,'vector');
UTrans = U';
permInv(perm) = 1:numPix;

x1k = b;
x2k = b;
uk = b;
y1k = applyB(x1k);
y2k = applyC(x2k);

w1k = zeros(size(uk));
w2k = zeros(size(uk));
z1k = zeros(size(y1k));
z2k = zeros(size(y2k));

costs = [];

figure('Name','uk inside ADMM')

for k = 1:maxIter
    
    % Minimize aug. Lagrangian w.r.t. (x1,x2).
    rhs1 = uk - w1k/t + applyBTrans(y1k - z1k/t);
    x1kp1 = ifft2(fft2(rhs1)./eigValsMtrx);
%     check = applyMtrx(x1kp1) - rhs1;
%     check = max(abs(check(:)));
    
    rhs2 = uk - w2k/t + applyCTrans(y2k - z2k/t);    
    rhs2 = rhs2(:); 
    temp = UTrans\rhs2(perm);    
    q = U\temp;    
    x2kp1 = q(permInv);
    x2kp1 = reshape(x2kp1,[numRows,numCols]);   
%     check = sparseMtrx*x2kp1(:) - rhs2;
%     check = max(abs(check(:)));
    
    % Minimize aug. Lagrangian w.r.t. (u,y1,y2);
    
    ukp1 = .5*(x1kp1 + w1k/t + x2kp1 + w2k/t);
    
    y1Hat = applyB(x1kp1) + z1k/t;
    y2Hat = applyC(x2kp1) + z2k/t;
    
    term = y1Hat + y2Hat;
    termData = term(:,:,1);
    termReg = term(:,:,2:end);
    
    tmpData = (b/beta^2 + t*termData/2)/(1/beta^2 + t/2);
%     tmpData = prox1Norm(termData - b,2/(t*beta)) + b;
%     tmpReg = proxIsoNorm(2*gamma/(t*beta),termReg);
    tmpReg = prox1Norm(termReg,2*gamma/(t*beta));
    tmp = cat(3,tmpData,tmpReg);    
    
    y1kp1 = (tmp + y1Hat - y2Hat)/2;
    y2kp1 = (tmp + y2Hat - y1Hat)/2;
    
%     w1kp1 = w1k + t*(x1kp1 - ukp1);
%     w2kp1 = w2k + t*(x2kp1 - ukp1);
%     z1kp1 = z1k + t*(applyB(x1kp1) - y1kp1);
%     z2kp1 = z2k + t*(applyC(x2kp1) - y2kp1);    
    
    w1kp1 = w1k + t*(overRelax*x1kp1 + (1-overRelax)*uk - ukp1);
    w2kp1 = w2k + t*(overRelax*x2kp1 + (1-overRelax)*uk - ukp1);
    z1kp1 = z1k + t*(overRelax*applyB(x1kp1) + (1-overRelax)*y1k - y1kp1);
    z2kp1 = z2k + t*(overRelax*applyC(x2kp1) + (1-overRelax)*y2k - y2kp1);     
    
    x1k = x1kp1;
    x2k = x2kp1;
    uk = ukp1;
    y1k = y1kp1;
    y2k = y2kp1;
    
    w1k = w1kp1;
    w2k = w2kp1;
    z1k = z1kp1;
    z2k = z2kp1;
    
    Dx = applyD(uk);
    KxMinusb = applyK(uk) - b;
    primalCost = (.5/beta^2)*sum(KxMinusb(:).^2) + (gamma/beta)*sum(abs(Dx(:)));
    primalCost = primalCost/gamma; % this is to compare with CP costs.
    
    costs = [costs,primalCost];         
                    
    if mod(k,showTrigger) == 0         
        imshow(uk,[])    
        disp(['admm iteration is: ',num2str(k)])
        disp(['primal cost is: ',num2str(primalCost)])         
        keyboard        
    end         
    
    
end

end















