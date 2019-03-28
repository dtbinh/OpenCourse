% This code can be used to reproduce the example in section 6.2,
% and in particular Figure 7, of the paper 
% "Image deblurring by primal-dual operator splitting."
% We solve minimize .5*||Kx - b||_F^2 + gamma*|| Dx ||_1.
% K = Kp + Ks, Ks is sparse, Kp uses periodic BCs.
% D is a tight frame analysis operator.

% maxIter = 2000;
maxIter = 300;
showTrigger = 10000; % Determines how frequently restored image is displayed.

testDR = 1; % primal-dual Douglas-Rachford
testCP = 1; % Chambolle-Pock
test_drPrimal = 1; % primal Douglas-Rachford
test_admm = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mtstream = RandStream('mt19937ar','seed',1);
RandStream.setGlobalStream(mtstream);

fname = 'barbara.png';
resizeFactor = .5;

I = imread(fname);
I = double(I(:,:,1));
I = imresize(I,resizeFactor);

[numRows,numCols] = size(I);
numPix = numRows*numCols;

mask = fspecial('gaussian',[9 9],4);
E = makeSparseA2_replicateBCs_deblur(numRows,numCols,mask);

sigmaNoise = (1e-3);
mu = 300000;
gamma = 1/mu;

mn = min(I(:));
I = I - mn;
mx = max(I(:));
I = I/mx;

figure('Name','image before blurring')
imshow(I,[])

noise = sigmaNoise*randn(numRows,numCols);
b = imfilter(I,mask,'replicate');

b = b + noise;
figure('Name','b')
imshow(b,[])


if test_admm == 1;         
                
    overRelax = 1.7;            
    t = .007; 
    beta = 1;              

    params_admm = struct('maxIter',maxIter,'t',t,'beta',beta,'showTrigger',showTrigger,...
                        'gamma',gamma,'overRelax',overRelax);                         

    tic                    
    [x_admm,costs_admm] = deblur_mixedTightFrame_admm(b,mask,E,params_admm);
    time_admm = toc;    
                        
end

if test_drPrimal == 1;     
                
    overRelax = 1.7;            
    t = 70; 
    beta = .8;               

    params_drPrimal = struct('maxIter',maxIter,'t',t,'beta',beta,'showTrigger',showTrigger,...
                        'gamma',gamma,'overRelax',overRelax);                         

    tic                    
    [x_drPrimal,costs_drPrimal] = deblur_mixedTightFrame_drPrimal(b,mask,E,params_drPrimal);
    time_drPrimal = toc;            
            
end

if testDR == 1;
                        
    overRelax = 1.5;    
    t = 1000;    
    s = 7/t;    

    paramsDR = struct('maxIter',maxIter,'t',t,'s',s,'showTrigger',showTrigger,...
                    'nu',1,'gamma',gamma,'overRelax',overRelax);                         

    tic             
    [xDR,costsDR] = deblur_mixedTightFrame_DR(b,mask,E,paramsDR);
    timeDR = toc;           
            
end

if testCP == 1;                        
        
    tau = .00008;                 
    overRelax = 1.7;

    paramsCP = struct('rho',mu,'maxIter',maxIter,'tau',tau,...
                'showTrigger',showTrigger,'overRelax',overRelax);                          

    tic                 
    [xCP,costsCP] = deblur_mixedTightFrame_CP_orVersion(b,mask,E,paramsCP);
    timeCP = toc;  
    
end
 