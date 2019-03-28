function [xPlus_x1,costs,funVals] = solveTestProb3_drPrimal(b,B,C,params)
% We solve minimize ||x||_2 + || (B + C)x - b ||_1. No constraints.

t = params.t; 
beta = params.beta;
showTrigger = params.showTrigger;
overRelax = params.overRelax;
maxIter = params.maxIter;
xTrue = params.xTrue;
weight = params.weight;

%%%
B = beta*B;
C = beta*C;
b = beta*b;

% Now we're minimizing || x ||_2 + (weight/beta)||(B + C)x - b||_1
%%%

[m,n] = size(B);

A = B + C;
mtrxB = B'*B + eye(n);
mtrxC = C'*C + eye(n);

RB = chol(mtrxB);
RBTrans = RB';
RC = chol(mtrxC);
RCTrans = RC';

p_x1 = zeros(n,1);
p_x2 = p_x1;
p_y1 = zeros(m,1);
p_y2 = zeros(m,1);

costs = [];
funVals = [];

for k = 1:maxIter
    
    xPlus_x1 = prox2Norm(.5*(p_x1 + p_x2),t/2);
    xPlus_x2 = xPlus_x1;
    
    z = prox1Norm(p_y1 + p_y2 - b,2*t/(beta/weight)) + b;
    xPlus_y1 = .5*(z + p_y1 - p_y2);
    xPlus_y2 = .5*(z + p_y2 - p_y1);
    
    x1Hat = 2*xPlus_x1 - p_x1;
    x2Hat = 2*xPlus_x2 - p_x2;
    y1Hat = 2*xPlus_y1 - p_y1;
    y2Hat = 2*xPlus_y2 - p_y2;
    
    rhsB = x1Hat + B'*y1Hat;
%     yPlus_x1Check = mtrxB\rhsB;
    yPlus_x1 = RB\(RBTrans\rhsB);
    yPlus_y1 = B*yPlus_x1;
    
    rhsC = x2Hat + C'*y2Hat;
%     yPlus_x2Check = mtrxC\rhsC;
    yPlus_x2 = RC\(RCTrans\rhsC);
    yPlus_y2 = C*yPlus_x2;
    
    pPlus_x1 = p_x1 + overRelax*(yPlus_x1 - xPlus_x1);
    pPlus_x2 = p_x2 + overRelax*(yPlus_x2 - xPlus_x2);
    pPlus_y1 = p_y1 + overRelax*(yPlus_y1 - xPlus_y1);
    pPlus_y2 = p_y2 + overRelax*(yPlus_y2 - xPlus_y2);
    
    p_x1 = pPlus_x1;
    p_x2 = pPlus_x2;
    p_y1 = pPlus_y1;
    p_y2 = pPlus_y2;  
        
    cost = norm(xPlus_x1 - xTrue);
    costs = [costs,cost];
    
    funVal = norm(xPlus_x1) + (weight/beta)*norm(A*xPlus_x1 - b,1);
    funVals = [funVals,funVal];
    
    if mod(k,showTrigger) == 0                  
        disp(['drPrimal iteration is: ',num2str(k)])
        disp(['primal cost is: ',num2str(cost)])         
        keyboard        
    end       
    
end



end