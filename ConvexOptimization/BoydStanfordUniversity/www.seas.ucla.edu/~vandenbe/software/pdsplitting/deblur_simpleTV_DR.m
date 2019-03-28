function [x,costs] = deblur_simpleTV_DR(b,kernel,params)

% This code corresponds to section 5.1 in the paper
% "Image deblurring by primal-dual operator splitting."                    
% We solve minimize mu*||Kx - b ||_1 + || Dx ||_{iso}
% subject to 0 <= x <= 1.
% K convolves with kernel.  D is a discrete gradient operator.

mu = params.mu; % The paper refers to gamma = 1/mu.
t = params.t; % t is the primal step size 
s = params.s; % s is the dual step size
overRelax = params.overRelax; % overRelax is called rho in the paper.
maxIter = params.maxIter;
showTrigger = params.showTrigger;

% t = tau;
beta = sqrt(s/t); % operator gets scaled by beta. (p. 7 in paper)

%%% Using two stepsizes is equivalent to modifying objective and the
%%% operators, and using a single step size t.
%%% K and D are multiplied by beta, and mu and b 
%%% are adjusted to compensate for that.
kernel = kernel*beta;
mu = mu/beta;
b = b*beta;
%%% Finished modifying objective.
% After this modification, we are now minimizing
% mu*||Kx - b||_1 + (1/beta)*||Dx||_{iso}  (using the notation of this code)
% We do this using Douglas-Rachford with a single step size t.

evalIsoNorm = @(y) sum(sum(sqrt(y(:,:,1).^2 + y(:,:,2).^2)));

[numRows,numCols] = size(b);

eigValArr_K = eigValArrForCyclicConvOp(kernel,numRows,numCols);
eigValArr_D1 = eigValArrForCyclicConvOp(beta*[-1 1]',numRows,numCols); % Note that D is being multiplied by beta
eigValArr_D2 = eigValArrForCyclicConvOp(beta*[-1 1],numRows,numCols); % as mentioned above.

eigValArr_KTrans = conj(eigValArr_K);
eigValArr_D1Trans = conj(eigValArr_D1);
eigValArr_D2Trans = conj(eigValArr_D2);

applyK = @(x) applyCyclicConv2D(x,eigValArr_K);
applyKTrans = @(x) applyCyclicConv2D(x,eigValArr_KTrans);
applyD1 = @(x) applyCyclicConv2D(x,eigValArr_D1);
applyD1Trans = @(x) applyCyclicConv2D(x,eigValArr_D1Trans);
applyD2 = @(x) applyCyclicConv2D(x,eigValArr_D2);
applyD2Trans = @(x) applyCyclicConv2D(x,eigValArr_D2Trans);

applyD = @(x) cat(3,applyD1(x),applyD2(x));
applyDTrans = @(y) applyD1Trans(y(:,:,1)) + applyD2Trans(y(:,:,2));

applyMtrx = @(x) (t^2)*applyKTrans(applyK(x)) + (t^2)*applyDTrans(applyD(x)) + x; 
eigValsMtrx = (t^2)*eigValArr_KTrans.*eigValArr_K ...
                + (t^2)*eigValArr_D1Trans.*eigValArr_D1 + (t^2)*eigValArr_D2Trans.*eigValArr_D2 + ones(numRows,numCols);                                                                                                            

p_x = b;
p_z1 = zeros(numRows,numCols);
p_z2 = zeros(numRows,numCols,2);

costs = [];
figure('Name','x inside Douglas-Rachford')

for k = 1:maxIter
                       
    % compute a resolvent at x, z1,z2.
    xPlus_x = min(p_x,1);
    xPlus_x = max(xPlus_x,0);    
    
    term = p_z1 - t*b;
    xPlus_z1 = min(term,mu);
    xPlus_z1 = max(xPlus_z1,-mu);    
 
    xPlus_z2 = projectOntoDualIsoNormUnitBall(beta*p_z2)/beta;     
                   
    term_x = 2*xPlus_x - p_x;
    term_z1 = 2*xPlus_z1 - p_z1;
    term_z2 = 2*xPlus_z2 - p_z2;
          
    rhs = term_x - t*applyKTrans(term_z1) - t*applyDTrans(term_z2);
    yPlus_x = ifft2(fft2(rhs)./eigValsMtrx);
    
%     check = applyMtrx(yPlus_x) - rhs;
%     check = max(abs(check(:)));
    
    yPlus_z1 = term_z1 + t*applyK(yPlus_x);
    yPlus_z2 = term_z2 + t*applyD(yPlus_x);       
    
    p_xPlus = p_x + overRelax*(yPlus_x - xPlus_x);
    p_z1Plus = p_z1 + overRelax*(yPlus_z1 - xPlus_z1);
    p_z2Plus = p_z2 + overRelax*(yPlus_z2 - xPlus_z2);
    
    p_x = p_xPlus;
    p_z1 = p_z1Plus;
    p_z2 = p_z2Plus;        
        
    KxMinusb = applyK(xPlus_x) - b;
    Dx = applyD(xPlus_x);
    cost = mu*sum(abs(KxMinusb(:))) + (1/beta)*evalIsoNorm(Dx);
    costs = [costs,cost];                
    
    if mod(k,showTrigger) == 0
        imshow(xPlus_x,[]) 
        disp(['Douglas-Rachford iteration is: ',num2str(k)])
        disp(['primal cost is: ',num2str(cost)])   
        keyboard        
    end             
                                           
end

x = xPlus_x;

end