function [x,costs] = deblur_simpleTightFrame_spingarn(b,kernel,params)

% This code corresponds to section 6.1 in the paper
% "Image deblurring by primal-dual operator splitting."
% We solve minimize mu*||Kx - b||_1 + || Dx ||_1
% subject to 0 <= x <= 1.
% K convolves with kernel using periodic boundary conditions.
% D is analysis operator for a shearlet tight frame.

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
% mu*||Kx - b||_1 + (1/beta)*||Dx||_1  subject to 0 <= x <= 1 (using the notation of this code)
% We do this using Douglas-Rachford with a single step size t.

[numRows,numCols] = size(b);

eigValArr_K = eigValArrForCyclicConvOp(kernel,numRows,numCols);
eigValArr_KTrans = conj(eigValArr_K);

applyK = @(x) applyCyclicConv2D(x,eigValArr_K);
applyKTrans = @(x) applyCyclicConv2D(x,eigValArr_KTrans);

scale = [3 3 3 4 4]; 
qmf1 = MakeONFilter('Symmlet',4);
qmf2 = MakeONFilter('Symmlet',4);
ndir = 0;
alpha = 2^(ndir+2)+2; % D^T D = alpha I.
alpha = (beta^2)*alpha; % This is because we have TWO STEPSIZES.
                        % When we modify D, alpha is modified also.

applyD = @(x) beta*shearletTransform(x,qmf1,qmf2,scale,ndir);
applyDTrans = @(y) beta*shearletAdjTransform(y,qmf1,qmf2,scale,ndir);

applyMtrx = @(x) applyKTrans(applyK(x)) + (alpha + 1)*x; 
eigValsMtrx = eigValArr_KTrans.*eigValArr_K + (alpha+1)*ones(numRows,numCols);                                                                                                            

p_x = b;
p_y1 = zeros(numRows,numCols);
p_y2 = applyD(p_y1);

costs = [];
figure('Name','x inside Spingarn')

for k = 1:maxIter
                       
    % compute a resolvent at p_x, p_y1,p_y2.
    xPlus_x = min(p_x,1);
    xPlus_x = max(xPlus_x,0);  
  
    xPlus_y1 = prox1Norm(p_y1 - b,mu*t) + b;
    xPlus_y2 = prox1Norm(p_y2,t/beta);        
                   
    term_x = 2*xPlus_x - p_x;
    term_y1 = 2*xPlus_y1 - p_y1;
    term_y2 = 2*xPlus_y2 - p_y2;
          
    rhs = term_x + applyKTrans(term_y1) + applyDTrans(term_y2);
    yPlus_x = ifft2(fft2(rhs)./eigValsMtrx);
    
%     check = applyMtrx(yPlus_x) - rhs;
%     check = max(abs(check(:)));

    yPlus_y1 = applyK(yPlus_x);
    yPlus_y2 = applyD(yPlus_x);          
    
    p_xPlus = p_x + overRelax*(yPlus_x - xPlus_x);
    p_y1Plus = p_y1 + overRelax*(yPlus_y1 - xPlus_y1);
    p_y2Plus = p_y2 + overRelax*(yPlus_y2 - xPlus_y2);
    
    p_x = p_xPlus;
    p_y1 = p_y1Plus;
    p_y2 = p_y2Plus;        
        
    KxMinusb = applyK(xPlus_x) - b;
    Dx = applyD(xPlus_x);
    cost = mu*sum(abs(KxMinusb(:))) + (1/beta)*sum(abs(Dx(:)));
    costs = [costs,cost];                
    
    if mod(k,showTrigger) == 0
        imshow(xPlus_x,[]) 
        disp(['Spingarn iteration is: ',num2str(k)])
        disp(['primal cost is: ',num2str(cost)])   
        keyboard        
    end             
                                           
end

x = xPlus_x;

end