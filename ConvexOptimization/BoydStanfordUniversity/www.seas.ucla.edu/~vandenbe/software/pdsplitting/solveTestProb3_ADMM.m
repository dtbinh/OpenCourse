function [uk,costs,funVals] = solveTestProb3_ADMM(b,B,C,params)

t = params.t; 
beta = params.beta;
showTrigger = params.showTrigger;
alpha = params.overRelax;
maxIter = params.maxIter;
xTrue = params.xTrue;
weight = params.weight;

%%%
B = beta*B;
C = beta*C;
b = beta*b;

% Now we're minimizing || x ||_1 + (weight/beta)||(B + C)x - b||_1
%%%

[m,n] = size(B);

A = B + C;
mtrxB = B'*B + eye(n);
mtrxC = C'*C + eye(n);

RB = chol(mtrxB);
RBTrans = RB';
RC = chol(mtrxC);
RCTrans = RC';

x1k = zeros(n,1);
x2k = x1k;
uk = x1k;
y1k = zeros(m,1);
y2k = zeros(m,1);

w1k = zeros(size(uk));
w2k = zeros(size(uk));
z1k = zeros(size(y1k));
z2k = zeros(size(y2k));

costs = [];
funVals = [];

for k = 1:maxIter
    
    rhsB = uk - w1k/t + B'*(y1k - z1k/t);
%     x1kp1Check = mtrxB\rhsB;
    x1kp1 = RB\(RBTrans\rhsB);
    
    rhsC = uk - w2k/t + C'*(y2k - z2k/t);
%     x2kp1Check = mtrxC\rhsC;
    x2kp1 = RC\(RCTrans\rhsC);
    
    ukp1 = prox2Norm(.5*(x1kp1 + w1k/t + x2kp1 + w2k/t),.5/t);
    
    y1Hat = B*x1kp1 + z1k/t;
    y2Hat = C*x2kp1 + z2k/t;
%     tmp = prox1Norm(y1Hat + y2Hat - b,2/(t*beta)) + b;
    tmp = prox1Norm(y1Hat + y2Hat - b,2/(t*beta/weight)) + b;
    y1kp1 = .5*(tmp + y1Hat - y2Hat);
    y2kp1 = .5*(tmp + y2Hat - y1Hat);
    
    w1kp1 = w1k + t*(alpha*x1kp1 + (1-alpha)*uk - ukp1);
    w2kp1 = w2k + t*(alpha*x2kp1 + (1-alpha)*uk - ukp1);
    z1kp1 = z1k + t*(alpha*(B*x1kp1) + (1-alpha)*y1k - y1kp1);
    z2kp1 = z2k + t*(alpha*(C*x2kp1) + (1-alpha)*y2k - y2kp1);     
    
    x1k = x1kp1;
    x2k = x2kp1;
    uk = ukp1;
    y1k = y1kp1;
    y2k = y2kp1;
    
    w1k = w1kp1;
    w2k = w2kp1;
    z1k = z1kp1;
    z2k = z2kp1;   
        
    cost = norm(xTrue - uk);
    costs = [costs,cost];
    
%     funVal = norm(uk,2) + (1/beta)*norm(A*uk - b,1);
    funVal = norm(uk,2) + (weight/beta)*norm(A*uk - b,1);
    funVals = [funVals,funVal];
    
    if mod(k,showTrigger) == 0                  
        disp(['ADMM iteration is: ',num2str(k)])
        disp(['cost is: ',num2str(cost)])         
        keyboard        
    end     
    
    
end



end








