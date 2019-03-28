% This code can be used to reproduce the example in section 7,
% and in particular Figure 10, of the paper 
% "Image deblurring by primal-dual operator splitting."
% We solve minimize || Kx - b||_1 + gamma*|| Dx ||_{iso}
% subject to  0 <= x <= 1
% K = Kbd + Ks. D = Dbd + Ds.
% D uses symmetric BCs.  K uses replicate BCs.
% Kbd and Dbd can be represented by block diagonal matrices 
% (with an appropriate ordering of pixels).

maxIter = 300;
% maxIter = 2000;
showTrigger = 10000; % Determines how frequently restored image is displayed.

testDR = 1; % primal-dual Douglas-Rachford
testCP = 1; % Chambolle-Pock
test_drPrimal = 1; % primal Douglas-Rachford
testADMM = 1;

nrMask = 9;
ncMask = 9;
R = 2;
S = 2;
kernels = zeros(nrMask,ncMask,R,S);
kernels(:,:,1,1) = fspecial('gaussian',[nrMask ncMask],4);
kernels(:,:,1,2) = fspecial('gaussian',[nrMask ncMask],3.5);
kernels(:,:,2,1) = fspecial('gaussian',[nrMask ncMask],3);
kernels(:,:,2,2) = fspecial('gaussian',[nrMask ncMask],2);

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
nr = numRows/R;
nc = numCols/S;

Ks = makeSparseA2_replicateBCs_varBlur(numRows,numCols,kernels);
Ks = Ks + makeSparseA2_varBlur(kernels,numRows,numCols);

eigValArrs_KpBlocks = zeros(nr,nc,R,S);
eigValArrs_KpTransBlocks = zeros(nr,nc,R,S);
for r = 1:R
    for s = 1:S
        kernel = kernels(:,:,r,s);
        arr = eigValArrForCyclicConvOp(kernel,nr,nc);
        eigValArrs_KpBlocks(:,:,r,s) = arr;
        eigValArrs_KpTransBlocks(:,:,r,s) = conj(arr);
    end
end

applyKp = @(x) applyBlockCyclicConv2D(x,eigValArrs_KpBlocks);
applyKpTrans = @(x) applyBlockCyclicConv2D(x,eigValArrs_KpTransBlocks);

applyKs = @(x) reshape(Ks*x(:),numRows,numCols);
applyKsTrans = @(x) reshape(Ks'*x(:),numRows,numCols);

applyK = @(x) applyKp(x) + applyKs(x);
applyKTrans = @(x) applyKpTrans(x) + applyKsTrans(x);

noiseDensity = .1;
gamma = 1/100;

mn = min(I(:));
I = I - mn;
mx = max(I(:));
I = I/mx;

figure('Name','image before blurring')
imshow(I,[])

b = applyK(I);

prm = randperm(numPix);
prm = prm(1:floor(noiseDensity*numPix));
vals = round(rand(size(prm)));
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
    [xADMM,costsADMM] = deblur_mixedTV_decomp_admm(b,kernels,Ks,paramsADMM);
    timeADMM = toc;
                        
end

if test_drPrimal == 1; 
                      
    overRelax = 1.9;            
    t = .4; 
    beta = 2;   

    params_drPrimal = struct('maxIter',maxIter,'t',t,'beta',beta,'showTrigger',showTrigger,...
                        'gamma',gamma,'overRelax',overRelax);                         

    tic                    
    [x_drPrimal,costs_drPrimal] = deblur_mixedTV_decomp_drPrimal(b,kernels,Ks,params_drPrimal);
    time_drPrimal = toc;
                        
end

if testDR == 1; 
                    
    overRelax = 1.9; 
    t = 1.2; 
    beta = 2.7;   

    paramsDR = struct('maxIter',maxIter,'t',t,'beta',beta,'showTrigger',showTrigger,...
                        'gamma',gamma,'overRelax',overRelax,'R',R,'S',S);                         

    tic                    
    [xDR,costsDR] = deblur_mixedTV_DR_varBlur(b,kernels,Ks,paramsDR);
    timeDR = toc;
                        
end

if testCP == 1; 
            
    tau = .03;
    overRelax = 1.9;

    paramsCP = struct('gamma',gamma,'tau',tau,'maxIter',maxIter,...
                    'showTrigger',showTrigger,'overRelax',overRelax);

    tic                 
    [xCP,costsCP] = deblur_mixedTV_decomp_CP_orVersion(b,kernels,Ks,paramsCP);   
    timeCP = toc; 
                              
end

