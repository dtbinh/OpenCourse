function [xPlus_x,costs,funVals] = solveTestProb3_mixedDR(b,B,C,params)
% We solve minimize ||x||_2 + || (B + C)x - b ||_1.  No constraints.

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

A = B + C;
[yDim,xDim] = size(B);

p_x = zeros(xDim,1);
p_y = zeros(yDim,1);
p_z = zeros(size(p_y));
p_w = zeros(size(p_x));

mtrxB = eye(xDim) + (t^2)*(B'*B);
mtrxC = (1 + t^2)*eye(xDim) + ((t^2)/(1 + t^2))*(C'*C);

RB = chol(mtrxB);
RBTrans = RB';
RC = chol(mtrxC);
RCTrans = RC';

costs = [];
funVals = [];

for k = 1:maxIter
   
    rhs1 = p_x - t*(B'*p_z);    
%     xPlus_x = mtrxB\rhs1;            
    xPlus_x = RB\(RBTrans\rhs1);
    
    xPlus_z = p_z + t*(B*xPlus_x);
        
    xPlus_y = prox1Norm(p_y - b,t/(beta/weight)) + b;
        
    prox = prox2Norm(p_w/t,1/t);
    xPlus_w = p_w - t*prox; % using Moreau decomposition.
    
    xHat = 2*xPlus_x - p_x;    
    yHat = 2*xPlus_y - p_y;
    zHat = 2*xPlus_z - p_z;
    wHat = 2*xPlus_w - p_w;  
    
    rhs2 = xHat + C'*(((t^2)/(1+t^2))*yHat - (t/(1+t^2))*zHat) - t*wHat;    
%     yPlus_xCheck = mtrxC\rhs2;
    yPlus_x = RC\(RCTrans\rhs2);
    
    CyPlus_x = C*yPlus_x;
    yPlus_y = (1/(1+t^2))*yHat + (t/(1+t^2))*zHat + ((t^2)/(1+t^2))*CyPlus_x;
    yPlus_z = (-t/(1+t^2))*yHat + (1/(1+t^2))*zHat + (t/(1+t^2))*CyPlus_x;
    yPlus_w = wHat + t*yPlus_x;    
    
    pPlus_x = p_x + overRelax*(yPlus_x - xPlus_x);
    pPlus_y = p_y + overRelax*(yPlus_y - xPlus_y);
    pPlus_z = p_z + overRelax*(yPlus_z - xPlus_z);
    pPlus_w = p_w + overRelax*(yPlus_w - xPlus_w);  
    
    p_x = pPlus_x;
    p_y = pPlus_y;
    p_z = pPlus_z;
    p_w = pPlus_w;            
    
    cost = norm(xPlus_x - xTrue);    
    costs = [costs,cost];       
    
    funVal = norm(xPlus_x) + (weight/beta)*norm(A*xPlus_x - b,1);
    funVals = [funVals,funVal];
                    
    if mod(k,showTrigger) == 0                    
        disp(['mixed splitting DR iteration is: ',num2str(k)])
        disp(['primal cost is: ',num2str(cost)])         
        keyboard        
    end      
    
end



end







