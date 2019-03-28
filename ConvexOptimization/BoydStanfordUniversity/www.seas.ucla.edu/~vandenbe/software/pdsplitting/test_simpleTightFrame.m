% This code can be used to reproduce the example in section 6.1,
% and in particular Figure 6, of the paper 
% We solve minimize mu*||Kx - b||_1 + || Dx ||_1
% subject to 0 <= x <= 1.
% K convolves with kernel using periodic boundary conditions.
% D is analysis operator for a shearlet tight frame.

maxIter = 300;
% maxIter = 2000;
showTrigger = 10000; % Determines how frequently restored image is displayed.

testCP = 1; % Chambolle-Pock
testADMM = 1;
testDR = 1; % primal-dual Douglas-Rachford
testSpingarn = 1; % primal Douglas-Rachford

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mtstream = RandStream('mt19937ar','seed',1);
RandStream.setGlobalStream(mtstream);

mask = fspecial('gaussian',[7 7],5);
noiseDensity = .3;
fname = 'cameraman.tif';

mu = 100; % gamma = 1/mu.
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

b = imfilter(I,mask,'circular');
b = imnoise(b,'salt & pepper',noiseDensity);
figure('Name','b')
imshow(b,[])


if testCP == 1;
    
    tau = .00005;
    overRelax = 1.7;

    paramsCP = struct('maxIter',maxIter,'tau',tau,'showTrigger',showTrigger,...
                    'rho',mu,'overRelax',overRelax); 
    tic
    [xCP,costsCP] = deblur_simpleTightFrame_CP_orVersion(b,mask,paramsCP);     
    timeCP = toc;
         
end

if testADMM == 1
                        
    stepSize = 400;
    beta = 1.5;
    alpha = 1.7;            

    paramsADMM = struct('rho',mu,'maxIter',maxIter,'stepSize',stepSize,...
                                    'showTrigger',showTrigger,'beta',beta,'gamma',alpha);

    tic    
    [xADMM,costsADMM] = deblur_simpleTightFrame_ADMM_overRelax(b,mask,paramsADMM);
    timeADMM = toc; 
             
end

if testDR == 1   
    
    overRelax = 1.9;    
    t = .8;    
    s = 1125;
    paramsDR = struct('mu',mu,'maxIter',maxIter,'t',t,'s',s,'showTrigger',showTrigger,'overRelax',overRelax);            
    
    tic            
    [xDR,costsDR] = deblur_simpleTightFrame_DR(b,mask,paramsDR);            
    timeDR = toc; 

end  

if testSpingarn == 1  
                                                    
    overRelax = 1.9;
    t = .0008;
    beta = 1;          
    s = (beta^2)*t;

    paramsSpingarn = struct('mu',mu,'t',t,'s',s,'overRelax',overRelax,'maxIter',maxIter,'showTrigger',showTrigger);

    tic
    [xSpingarn,costsSpingarn] = deblur_simpleTightFrame_spingarn(b,mask,paramsSpingarn);                
    timeSpingarn = toc;    
              
end 



