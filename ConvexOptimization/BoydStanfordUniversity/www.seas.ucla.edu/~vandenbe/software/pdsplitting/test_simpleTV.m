% This code can be used to reproduce the example in section 5.1,
% and in particular Figure 2, of the paper 
% "Image deblurring by primal-dual operator splitting."
% We solve minimize mu*||Kx - b ||_1 + || Dx ||_{iso}
% subject to 0 <= x <= 1.
% K convolves with kernel (see below).  D is a discrete gradient operator.

maxIter = 100;
showTrigger = 10000; % Determines how frequently restored image is displayed.

testADMM = 1;
testDR = 1; % primal-dual Douglas-Rachford
testCP = 1; % Chambolle-Pock
testSpingarn = 1; % primal Douglas-Rachford

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mtstream = RandStream('mt19937ar','seed',1);
RandStream.setGlobalStream(mtstream);

kernel = fspecial('gaussian',[15 15],7);
noiseDensity = .5;
fname = 'manWithHat.tiff';

mu = 10; % The paper refers to gamma = 1/mu
resizeFactor = 1;

I = imread(fname);
I = double(I(:,:,1));
I = imresize(I,resizeFactor);

[numRows,numCols] = size(I);
numPix = numRows*numCols;

mn = min(I(:));
I = I - mn;
mx = max(I(:));
I = I/mx;

figure('Name','image before blurring')
imshow(I,[])

b = imfilter(I,kernel); % Note that zero boundary conditions are being used here, 
                      % whereas periodic BCs are used in reconstruction.  
b = imnoise(b,'salt & pepper',noiseDensity);
figure('Name','b')
imshow(b,[])

if testCP == 1;  
            
    overRelax = 1.7;
    tau = .003;

    paramsCP = struct('rho',mu,'maxIter',maxIter,'tau',tau,...
                'showTrigger',showTrigger,'overRelax',overRelax);

    tic
    [xCP,costsCP] = deblur_simpleTV_CP_orVersion(b,kernel,paramsCP);       
    timeCP = toc;          
        
end

if testADMM== 1    

    alpha = 1.8; %overrelaxation parameter.
    beta = 1; % problem scaling parameter.
    stepSize = 50;
    paramsADMM = struct('rho',mu,'maxIter',maxIter,'stepSize',stepSize,...
                                    'showTrigger',showTrigger,'beta',beta,'alpha',alpha);

    tic    
    [xADMM,costsADMM] = deblur_simpleTV_ADMM_overRelax(b,kernel,paramsADMM);
    timeADMM = toc; 
             
end

if testDR == 1          
                              
    overRelax = 1.9;
    t = 1.5;
    s = 150;

    paramsDR = struct('mu',mu,'t',t,'s',s,'overRelax',overRelax,'maxIter',maxIter,'showTrigger',showTrigger);

    tic
    [xDR,costsDR] = deblur_simpleTV_DR(b,kernel,paramsDR);                
    timeDR = toc;
              
end  


if testSpingarn == 1  
                        
    overRelax = 1.9;
    t = .008;
    beta = 1.2;
    s = (beta^2)*t;

    paramsSpingarn = struct('mu',mu,'t',t,'s',s,'overRelax',overRelax,'maxIter',maxIter,'showTrigger',showTrigger);

    tic
    [xSpingarn,costsSpingarn] = deblur_simpleTV_spingarn(b,kernel,paramsSpingarn);                
    timeSpingarn = toc;    
              
end 

