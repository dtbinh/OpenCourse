function [xPlus_x,costs] = deblur_mixedTightFrame_DR(b,kernel,Ks,params)

% This code corresponds to section 6.2 in the paper
% "Image deblurring by primal-dual operator splitting." 
% We solve minimize .5*nu*||Kx - b||_F^2 + gamma*|| Dx ||_1.
% K = Kp + Ks, Ks is sparse, Kp uses periodic BCs.
% D is a tight frame analysis operator.

maxIter = params.maxIter;
t = params.t; % primal step size
s = params.s; % dual step size
nu = params.nu;
gamma = params.gamma;
showTrigger = params.showTrigger;
overRelax = params.overRelax; % called rho in the paper.

gamma_orig = gamma;

beta = sqrt(s/t); % operator gets scaled by beta. (p. 7 in paper)

%%% Using two stepsizes is equivalent to modifying objective and the
%%% operators.  See p. 7 in paper.
Ks = beta*Ks;
kernel = beta*kernel;
nu = nu/beta^2;
gamma = gamma/beta;
b = beta*b;
%%% Finished modifying objective.

[numRows,numCols] = size(b);
numPix = numRows*numCols;

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

factor1 = (t^2)/(1 + t/nu);
factor2 = 1 + (t^2)*alpha;
applyMtrx = @(x) factor1*applyKpTrans(applyKp(x)) + factor2*x; 
eigValsMtrx = factor1*eigValArr_KpTrans.*eigValArr_Kp ...
                + factor2*ones(numRows,numCols);

sparseMtrx = [speye(numPix),t*Ks';-t*Ks,speye(numPix)];
[L,U,P_perm,Q_perm] = lu(sparseMtrx);

p_x = b;
p_w = zeros(size(b));
p_z = applyD(p_w);

costs = [];
figure('Name','xPlus_x inside DR')

for n = 1:maxIter      
        
    rhs = p_x - t*applyKpTrans((p_w - t*b)/(1+t/nu)) - t*applyDTrans(p_z);
    xPlus_x = ifft2(fft2(rhs)./eigValsMtrx);            
    
    xPlus_w = (p_w - t*b + t*applyKp(xPlus_x))/(1+t/nu);
    xPlus_z = p_z + t*applyD(xPlus_x);
    
    xHat = 2*xPlus_x - p_x;
    wHat = 2*xPlus_w - p_w;
    zHat = 2*xPlus_z - p_z;
    
    rhs2 = cat(3,xHat,wHat);
    rhs2 = rhs2(:);
    
    Prhs2 = P_perm*rhs2;
    temp = L\Prhs2;
    w = U\temp;
    both = Q_perm*w;
    both = reshape(both,[numRows,numCols,2]);  
    
    yPlus_x = both(:,:,1);
    yPlus_w = both(:,:,2);    
    
    yPlus_z = min(zHat,gamma);
    yPlus_z = max(yPlus_z,-gamma);
    
    pPlus_x = p_x + overRelax*(yPlus_x - xPlus_x);
    pPlus_w = p_w + overRelax*(yPlus_w - xPlus_w);
    pPlus_z = p_z + overRelax*(yPlus_z - xPlus_z);
    
    p_x = pPlus_x;
    p_w = pPlus_w;
    p_z = pPlus_z;
    
    Dx = applyD(xPlus_x);
    KxMinusb = applyK(xPlus_x) - b;
    primalCost = (.5*nu)*sum(KxMinusb(:).^2) + gamma*sum(abs(Dx(:)));
    primalCost = primalCost/gamma_orig; % this is to compare with CP costs.
    
    costs = [costs,primalCost];      
                    
    if mod(n,showTrigger) == 0        
        imshow(xPlus_x,[])  
        disp(['Douglas-Rachford iteration is: ',num2str(n)])
        disp(['primal cost is: ',num2str(primalCost)])         
        keyboard         
    end             
    
end


end