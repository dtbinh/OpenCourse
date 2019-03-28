function [xPlus_x,costs,errors] = tvDeblur_varBlur_freeBCs_DR_huber(b,kernels,U,params,xStar)

% b has size M by N.
% kernels have size (1 + 2d) by (1 + 2d).
% bHat has size (M + 2d) by (N + 2d).
% x will have size (M + 2d) by (N + 2d).
% E simply zeros out border pixels, leaving other pixels unchanged.
% We solve minimize huber(Kx - b) + gamma*||Dx||_{iso}
% K = U1*K1 + ... + U_P*K_P.
% Each Kp convolves using periodic BCs.
% D is a discreet gradient operator using periodic boundary conditions.
% Ui zeros out pixels on the border, and multiplies remaining components
% according to how much Ki should contribute to the output of K at that
% component.
% Ui will be stored in the code as an array of size M by N, rather than an array of
% size M + 2d by N + 2d, because there's no need to store a bunch of zeros.

maxIter = params.maxIter;
showTrigger = params.showTrigger;
t = params.t; % primal step size
beta = params.beta;
gamma = params.gamma;
overRelax = params.overRelax; 
mu = params.mu;  % Huber penalty parameter.

[numRows,numCols,numKernels] = size(U);
d = (size(kernels,1) - 1)/2 ; % We assume the kernels have size (1 + 2d) by (1 + 2d).

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Account for beta. %%%%
kernels = beta*kernels;
b = beta*b;
gamma = gamma/beta;

eigValArr_D1 = eigValArrForCyclicConvOp(beta*[-1 1]',numRows + 2*d,numCols + 2*d); % Note that D is being multiplied by beta
eigValArr_D2 = eigValArrForCyclicConvOp(beta*[-1 1],numRows + 2*d,numCols + 2*d); % as mentioned above.

%%% Now we're minimizing .5|| (Kx - b)/beta||^2 + gamma ||Dx||_1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eigValArr_D1Trans = conj(eigValArr_D1);
eigValArr_D2Trans = conj(eigValArr_D2);

applyD1 = @(x) applyCyclicConv2D(x,eigValArr_D1);
applyD1Trans = @(x) applyCyclicConv2D(x,eigValArr_D1Trans);
applyD2 = @(x) applyCyclicConv2D(x,eigValArr_D2);
applyD2Trans = @(x) applyCyclicConv2D(x,eigValArr_D2Trans);

applyD = @(x) cat(3,applyD1(x),applyD2(x));
applyDTrans = @(y) applyD1Trans(y(:,:,1)) + applyD2Trans(y(:,:,2));

evalIsoNorm = @(y) sum(sum(sqrt(y(:,:,1).^2 + y(:,:,2).^2)));

Lambda = (sum(U.^2,3)).^(-1);
v = zeros(numRows,numCols,numKernels);
for k = 1:numKernels
    v(:,:,k) = U(:,:,k).*Lambda;
end
vNrmSquared = sum(v.^2,3);

eigValArrs_Kp = zeros(numRows + 2*d,numCols + 2*d,numKernels);
eigValArrs_KpTrans = zeros(numRows + 2*d,numCols + 2*d,numKernels);
for k = 1:numKernels
    eigValArr = eigValArrForCyclicConvOp(kernels(:,:,k),numRows + 2*d,numCols + 2*d);
    eigValArrs_Kp(:,:,k) = eigValArr;
    eigValArrs_KpTrans(:,:,k) = conj(eigValArr);
end

% mtrx = I + (t^2)A^T A.
eigValsMtrx = zeros(numRows + 2*d,numCols + 2*d);
for k = 1:numKernels
    eigValsMtrx = eigValsMtrx + eigValArrs_KpTrans(:,:,k).*eigValArrs_Kp(:,:,k);
end
eigValsMtrx = eigValsMtrx + eigValArr_D1Trans.*eigValArr_D1 + eigValArr_D2Trans.*eigValArr_D2;
eigValsMtrx = (t^2)*eigValsMtrx + ones(numRows + 2*d,numCols + 2*d);

z_x = zeros(numRows + 2*d,numCols + 2*d);
% z_x(d+1:d + numRows,d+1:d + numCols) = b;
z_x(d+1:d + numRows,d+1:d + numCols) = b/beta;

z_z1 = zeros(numRows + 2*d,numCols + 2*d,numKernels);
z_z2 = zeros(size(applyD(z_x)));

costs = [];
errors = [];

for iter = 1:maxIter
    
    xPlus_x = z_x; % prox operator of zero function is identity.
    
    % Evaluate prox operator of g1Star at z_z1, using Moreau decomposition.
    % First evaluate proxg1 at z_z1/t , with step size 1/t.
    yHat = z_z1/t; 
    
    %%%%%%%%%%%%%%%%%%%%%%% Compute proxg1 using eqs 15 and 16 in paper. %%%%%%%%%%%%%%%%%

    UyHat = zeros(numRows,numCols);
    for k = 1:numKernels
        UyHat = UyHat + U(:,:,k).*yHat(d+1:d+numRows,d+1:d+numCols,k);
    end    

    in = UyHat - b;
    stepSize = (1/t)*(Lambda.^(-1));    
    
%     proxSpecial = in./(stepSize/beta^2 + 1);    % prox of function phi_f^i(x_i/beta) evaluated at in.     
    proxSpecial = beta*proxHuber(in/beta,mu,stepSize/beta^2);
    
    foo = Lambda.*(proxSpecial + b - UyHat);
    proxg1 = zeros(numRows + 2*d,numCols + 2*d,numKernels); 
    for k = 1:numKernels
        proxg1(d+1:d+numRows,d+1:d+numCols,k) = U(:,:,k).*foo;
    end 
    proxg1 = yHat + proxg1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % This was the original code to compute proxg1.

%     UyHat = zeros(numRows,numCols);
%     for k = 1:numKernels
%         UyHat = UyHat + U(:,:,k).*yHat(d+1:d+numRows,d+1:d+numCols,k);
%     end
%     LambdaUyHat = Lambda.*UyHat;
%     w = zeros(numRows,numCols,numKernels);
%     for k = 1:numKernels
%         w(:,:,k) = U(:,:,k).*LambdaUyHat;
%     end
%         
%     term = sum(v.*w,3)./vNrmSquared;    
%     
%     tau = (t*vNrmSquared).^(-1); 
%     tmp = (tau.*b + term)./(tau + 1);    
%     vec = Lambda.*(tmp - UyHat);
%     proxg1 = zeros(numRows + 2*d,numCols + 2*d,numKernels);    
%     for k = 1:numKernels
%         proxg1(d+1:d+numRows,d+1:d+numCols,k) = U(:,:,k).*vec;
%     end
%     proxg1 = yHat + proxg1;    % Note that border components are handled correctly here.
    %%%%%%%%
    
    proxg1Star = z_z1 - t*proxg1;
    
%     proxg1StarCheck = checkProxg1Star(z_z1,t,U,b);
    
%     proxg1StarCheck = computeProxg1Star_ex1(z_z1,t,U,b);
%     diff = proxg1Star - proxg1StarCheck;
%     diff = max(abs(diff(:)));
    
%     proxg2Star = min(z_z2,gamma);
%     proxg2Star = max(proxg2Star,-gamma);
    
    proxg2Star = gamma*projectOntoDualIsoNormUnitBall(z_z2/gamma);    
    
    xPlus_z1 = proxg1Star;
    xPlus_z2 = proxg2Star;
    
    term_x = 2*xPlus_x - z_x;
    term_z1 = 2*xPlus_z1 - z_z1;
    term_z2 = 2*xPlus_z2 - z_z2;
    
    ATrans_term_z = zeros(numRows + 2*d,numCols + 2*d);
    for k = 1:numKernels
        ATrans_term_z = ATrans_term_z + applyCyclicConv2D(term_z1(:,:,k),eigValArrs_KpTrans(:,:,k));
    end    
    
    ATrans_term_z = ATrans_term_z + applyDTrans(term_z2);
    
    rhs = term_x - t*ATrans_term_z;
    yPlus_x = ifft2(fft2(rhs)./eigValsMtrx); 
    
    yPlus_z1 = zeros(numRows + 2*d,numCols + 2*d,numKernels);
    for k = 1:numKernels
        yPlus_z1(:,:,k) = applyCyclicConv2D(yPlus_x,eigValArrs_Kp(:,:,k));
    end    
    
    yPlus_z1 = term_z1 + t*yPlus_z1;
    yPlus_z2 = term_z2 + t*applyD(yPlus_x);
    
    z_xPlus = z_x + overRelax*(yPlus_x - xPlus_x);
    z_z1Plus = z_z1 + overRelax*(yPlus_z1 - xPlus_z1);
    z_z2Plus = z_z2 + overRelax*(yPlus_z2 - xPlus_z2);
    
    z_x = z_xPlus;
    z_z1 = z_z1Plus;
    z_z2 = z_z2Plus;      
    
               
    Kx = zeros(numRows,numCols);
    for k = 1:numKernels

        Kpx = applyCyclicConv2D(xPlus_x,eigValArrs_Kp(:,:,k));
        Kx = Kx + U(:,:,k).*Kpx(d+1:d+numRows,d+1:d+numCols);

    end
    Dx = applyD(xPlus_x);
    KxMinusb = Kx - b;
    prx = prox1Norm(KxMinusb/beta, mu);
    huberTerm = abs(prx) + (.5/mu)*(prx - KxMinusb/beta).^2;
    huberTerm = sum(huberTerm(:));
%     cost = (.5/beta^2)*sum(sum((Kx - b).^2)) + gamma*evalIsoNorm(Dx);
    cost = huberTerm + gamma*evalIsoNorm(Dx);
%     cost = cost/beta^2;
    costs = [costs,cost];

    if nargin == 5
        error = xPlus_x - xStar;
        error = norm(error(:));        
        errors = [errors,error];
    end
    
    if mod(iter,showTrigger) == 0
        
        if iter == showTrigger
            figure('Name','x inside primal-dual DR')
        end
        
        reduced = xPlus_x(d + 1:d + numRows, d + 1:d + numCols);
        
        imshow(reduced,[]) 
        disp(['primal-dual DR iteration is: ',num2str(iter)])
        disp(['primal cost is: ',num2str(cost)])   
%         keyboard     
        
    end     
    
end

% reduced = xPlus_x(d + 1:d + numRows, d + 1:d + numCols);


end