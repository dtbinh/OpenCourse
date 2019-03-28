% This code can be used to reproduce the example in section 5.2,
% and in particular Figure 4, of the paper 
% "Image deblurring by primal-dual operator splitting."
% We solve minimize || Kx - b||_1 + gamma*|| Dx ||_{iso}
% subject to  0 <= x <= 1
% K = Kp + Ks. D = Dp + Ds.
% D uses symmetric BCs.

maxIter = 300;
% maxIter = 2000;
showTrigger = 10000; % Determines how frequently restored image is displayed.

testDR = 1; % primal-dual Douglas-Rachford
testCP = 1; % Chambolle-Pock
test_drPrimal = 1; % primal Douglas-Rachford
testADMM = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mtstream = RandStream('mt19937ar','seed',1);
RandStream.setGlobalStream(mtstream);

fname = 'cameraman.tif';
resizeFactor = 1;
I = imread(fname);
I = double(I(:,:,1));
I = imresize(I,resizeFactor);

[numRows,numCols] = size(I);
numPix = numRows*numCols;

mask = fspecial('gaussian',[9 9],4);
E = makeSparseA2_replicateBCs_deblur(numRows,numCols,mask);

noiseDensity = .1;
gamma = 1/100;

mn = min(I(:));
I = I - mn;
mx = max(I(:));
I = I/mx;

figure('Name','image before blurring')
imshow(I,[])

b = imfilter(I,mask,'replicate');

prm = randperm(numPix);
prm = prm(1:floor(noiseDensity*numPix));
vals = sum(mask(:))*round(rand(size(prm)));
b(prm) = vals;
figure('Name','b')
imshow(b,[])

if testADMM == 1; 
                     
    overRelax = 1.6;  
    t = 1.5; 
    beta = 1.7;   

    paramsADMM = struct('maxIter',maxIter,'t',t,'beta',beta,'showTrigger',showTrigger,...
                        'gamma',gamma,'overRelax',overRelax);                         

    tic                    
    [xADMM,costsADMM] = deblur_mixedTV_admm(b,mask,E,paramsADMM);
    timeADMM = toc;    
                        
end

if test_drPrimal == 1; 
                      
    overRelax = 1.9;
    t = .4; 
    beta = 2;   

    params_drPrimal = struct('maxIter',maxIter,'t',t,'beta',beta,'showTrigger',showTrigger,...
                        'gamma',gamma,'overRelax',overRelax);                         

    tic                    
    [x_drPrimal,costs_drPrimal] = deblur_mixedTV_drPrimal(b,mask,E,params_drPrimal);
    time_drPrimal = toc;
                                    
end


if testDR == 1;       
    
    overRelax = 1.9;      
    t = 1.5;    
    s = 25/t;

    paramsDR = struct('maxIter',maxIter,'t',t,'s',s,'showTrigger',showTrigger,...
                        'gamma',gamma,'overRelax',overRelax);                         

    tic                    
    [xDR,costsDR] = deblur_mixedTV_DR(b,mask,E,paramsDR);
    timeDR = toc;
            
end


if testCP == 1;              
        
    tau = .02;
    overRelax = 1.9;

    paramsCP = struct('gamma',gamma,'tau',tau,'maxIter',maxIter,...
                'showTrigger',showTrigger,'theta',1,'overRelax',overRelax);

    tic                 
    [xCP,costsCP] = deblur_mixedTV_CP_orVersion(b,mask,E,paramsCP);         
    timeCP = toc;                                
                      
end

