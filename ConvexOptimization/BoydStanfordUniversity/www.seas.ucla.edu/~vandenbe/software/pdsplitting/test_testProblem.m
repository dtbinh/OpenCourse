% This code can be used to reproduce the example in section 3.5,
% and in particular Figure 1, of the paper 
% "Image deblurring by primal-dual operator splitting."
% We solve minimize ||x||_2 + ||(B + C)x - b||_1 

rng(9); % tuned for this.

showTrigger = 100000;
maxIter = 600;
weight = 1/100;

xDim = 500;
yDim = 500;

B = randn(yDim,xDim);
C = randn(yDim,xDim);
b = randn(yDim,1);

A = B + C;

test_ADMM = 1;
test_mixedDR = 1; % primal-dual Douglas-Rachford
test_drPrimal = 1; % primal Douglas-Rachford
test_cvx = 1;

if test_cvx == 1
        
    cvx_begin
        cvx_precision best
        variable x(xDim)
        minimize norm(x) + weight*norm(A*x - b,1)
    cvx_end
    xCvx = x;
    
    fStar = cvx_optval;
            
end

if test_ADMM == 1
                                  
    overRelax = 1.9;
    t = 1.7;
    beta = .05;                        

    params_ADMM = struct('t',t,'beta',beta,'showTrigger',showTrigger','maxIter',...
                                maxIter,'overRelax',overRelax,'xTrue',xCvx,'weight',weight);

    tic
    [x_ADMM,costs_ADMM,fVals_ADMM] = solveTestProb3_ADMM(b,B,C,params_ADMM);
    timeADMM = toc;
                     
end

if test_drPrimal == 1
            
    overRelax = 1.8; 
    t = .85;
    beta = .1;

    params_drPrimal = struct('t',t,'beta',beta,'showTrigger',showTrigger','maxIter',...
                                maxIter,'overRelax',overRelax,'xTrue',xCvx,'weight',weight);

    tic
    [x_drPrimal,costs_drPrimal,fVals_drPrimal] = solveTestProb3_drPrimal(b,B,C,params_drPrimal);
    time_drPrimal = toc;
                    
end

if test_mixedDR == 1
                                    
    overRelax = 1.9;  
    t = .8; 
    beta = .1;


    params_mixedDR = struct('t',t,'beta',beta,'showTrigger',showTrigger','maxIter',...
                                maxIter,'overRelax',overRelax,'xTrue',xCvx,'weight',weight);

    tic
    [x_mixedDR,costs_mixedDR,fVals_mixedDR] = solveTestProb3_mixedDR(b,B,C,params_mixedDR);
    time_mixedDR = toc;
                
end

figure
semilogy(costs_ADMM)
hold all
semilogy(costs_drPrimal)
semilogy(costs_mixedDR)
legend('ADMM','primal DR','primal-dual DR')


disp('finished')