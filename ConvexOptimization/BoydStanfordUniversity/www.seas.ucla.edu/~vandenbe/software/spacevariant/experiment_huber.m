
% This .m file reproduces the result of section 5.1 in 
% the paper "Total variation image deblurring
% with space-varying kernel", D. O'Connor and L. Vandenberghe.

rng(1);

testCP = 1;
testPrimalDualDR = 1;
maxIter = 1000;
maxIter = 500;

showTrigger = 100000;
load 'xDR_10k_huber.mat' % Restored image after 10k iters of DR.
load 'minPrimal_barb_huber_DR10k.mat' % Obj. fun value after 10k iters of DR.

R = 2; 
S = 2;
numKernels = R*S;

% Create kernels.
stdDevs = (1:numKernels);
nrMask = 17;
ncMask = nrMask;

kernels = zeros(nrMask,ncMask,numKernels);
for k = 1:numKernels
    kernels(:,:,k) = fspecial('gaussian',[nrMask,ncMask],stdDevs(k));
end

I = imread('barbara.png');
I = double(I(:,:,1));

d = (nrMask - 1)/2;

[numRows,numCols] = size(I);
numPix = numRows*numCols;

sigmaNoise = (1e-3);

gamma = .05;
mu = .001;

mn = min(I(:));
I = I - mn;
mx = max(I(:));
I = I/mx;

figure('Name','image before blurring')
imshow(I,[])

% Create U.
nr = numRows/R;
nc = numCols/S;
U = zeros(numRows,numCols,numKernels);
h = fspecial('gaussian',[17 17],8);
k = 1;
quads = [2 1 3 4];
for r = 1:R
    for s = 1:S
        
        rowStart = nr*(r-1) + 1;
        rowStop = nr*(r-1) + nr;
        colStart = nc*(s-1) + 1;
        colStop = nc*(s-1) + nc;
        
        arr = zeros(numRows,numCols);
        arr(rowStart:rowStop,colStart:colStop) = 1;
        arr = imfilter(arr,h,'replicate');
        U(:,:,quads(k)) = arr;
        k = k + 1;
        
    end
end
total = sum(U,3);
for k = 1:numKernels
    U(:,:,k) = U(:,:,k)./total;
end

% Create the blurry image b.
b = zeros(numRows,numCols);
for k = 1:numKernels
    im = imfilter(I,kernels(:,:,k),'replicate');
    b = b + U(:,:,k).*im;
end 

noise = sigmaNoise*randn(numRows,numCols);
b = b + noise;

noiseDensity = .1;
b = imnoise(b,'salt & pepper',noiseDensity);

figure('Name','b')
imshow(b,[])

if testCP == 1;             

    tau = .09;
    overRelax = 1.8;

    paramsCP = struct('gamma',gamma,'maxIter',maxIter,'tau',tau,...
                        'showTrigger',showTrigger,'overRelax',overRelax,'mu',mu);

    tic                 
    [xCP,costsCP,errorsCP] = tvDeblur_varBlur_freeBCs_CP_huber(b,kernels,U,paramsCP,xStar);
    timeCP = toc;               
                
end


if testPrimalDualDR == 1
     
    t = 900;
    beta = .05;
    overRelax = 1.9;

    paramsPDDR = struct('maxIter',maxIter,'showTrigger',showTrigger,'t',t,'beta',beta,...
                        'gamma',gamma,'overRelax',overRelax,'mu',mu);    
    tic
    [xPDDR,costsPDDR,errorsPDDR] = tvDeblur_varBlur_freeBCs_DR_huber(b,kernels,U,paramsPDDR,xStar);
    timePDDR = toc;     
                  
end

if testCP == 1 & testPrimalDualDR == 1

    fontSize = 18;    
    
    figure
    semilogy((costsCP - minPrimal)/minPrimal)
    hold on
    semilogy((costsPDDR - minPrimal)/minPrimal,'color','red')  
    title('$(F^k - F^\star)/F^\star$','Interpreter','Latex','FontSize',fontSize)
    xlabel('iteration number $k$','Interpreter','Latex','FontSize',fontSize)   
    legend('CP','DR')
   

    nrmxStar = norm(xStar(:));
    figure
    semilogy(errorsCP/nrmxStar)
    hold on
    semilogy(errorsPDDR/nrmxStar,'color','red')
    title('$\|x^k - x^\star \|/\|x^\star \|$','Interpreter','Latex','FontSize',fontSize)
    xlabel('iteration number $k$','Interpreter','Latex','FontSize',fontSize)   
    legend('CP','DR')
    
end

