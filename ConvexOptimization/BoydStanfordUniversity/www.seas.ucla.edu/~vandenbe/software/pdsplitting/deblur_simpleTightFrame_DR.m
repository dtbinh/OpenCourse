function [x,costs] = deblur_simpleTightFrame_DR(b,kernel,params)

% This code corresponds to section 6.1 in the paper
% "Image deblurring by primal-dual operator splitting."
% We solve minimize mu*||Kx - b||_1 + || Dx ||_1
% subject to 0 <= x <= 1.
% K convolves with kernel using periodic boundary conditions.
% D is analysis operator for a shearlet tight frame.

maxIter = params.maxIter;
showTrigger = params.showTrigger;
t = params.t; % primal step size
s = params.s; % dual step size
overRelax = params.overRelax; % called rho in paper
mu = params.mu; % The paper refers to gamma = 1/mu.

%%%%%%%
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

applyD = @(x) shearletTransform(x,qmf1,qmf2,scale,ndir);
applyDTrans = @(y) shearletAdjTransform(y,qmf1,qmf2,scale,ndir);

applyMtrx = @(x) (s*t)*applyKTrans(applyK(x)) + (s*t*alpha + 1)*x; 
eigValsMtrx = (s*t)*eigValArr_KTrans.*eigValArr_K ...
                + (s*t*alpha + 1)*ones(numRows,numCols);

%%%%%%%%%%                                                                                                                   

p_x = b;
p_z1 = zeros(numRows,numCols);
p_z2 = zeros(numRows,numCols,alpha);

costs = [];
figure('Name','xPlus_x inside DR')

for k = 1:maxIter
                       
    % compute a resolvent at p_x, p_z1, p_z2.
       
    xPlus_x = min(p_x,1);
    xPlus_x = max(xPlus_x,0);
        
    term = p_z1 - s*b;
    xPlus_z1 = min(term,mu);
    xPlus_z1 = max(xPlus_z1,-mu);         
      
    xPlus_z2 = min(p_z2,1);
    xPlus_z2 = max(xPlus_z2,-1);                
               
    term_x = 2*xPlus_x - p_x;
    term_z1 = 2*xPlus_z1 - p_z1;
    term_z2 = 2*xPlus_z2 - p_z2;
    
    rhs = term_x - t*applyKTrans(term_z1) - t*applyDTrans(term_z2);
    yPlus_x = ifft2(fft2(rhs)./eigValsMtrx);
    
%     check = applyMtrx(yPlus_x) - rhs;
%     check = max(abs(check(:)));
    
    yPlus_z1 = term_z1 + s*applyK(yPlus_x);
    yPlus_z2 = term_z2 + s*applyD(yPlus_x);            
    
    p_xPlus = p_x + overRelax*(yPlus_x - xPlus_x);
    p_z1Plus = p_z1 + overRelax*(yPlus_z1 - xPlus_z1);
    p_z2Plus = p_z2 + overRelax*(yPlus_z2 - xPlus_z2);
    
    p_x = p_xPlus;
    p_z1 = p_z1Plus;
    p_z2 = p_z2Plus;    
    
    KxPlus_xMinusb = applyK(xPlus_x) - b;
    DxPlus_x = applyD(xPlus_x);     
    cost = mu*sum(abs(KxPlus_xMinusb(:))) + sum(abs(DxPlus_x(:)));            
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