
% This .m file reproduces the result of section 5.2 in 
% the paper "Total variation image deblurring
% with space-varying kernel", D. O'Connor and L. Vandenberghe.

% This code will take about three minutes to run.

maxIter = 250;
showTrigger = 500;

load 'motorcycleSeg.mat'
figure('Name','Image with motion blur')
imshow(img,[])
figure('Name','segmentation')
imshow(lblimg,[])

% The segmentation stored in the variable lbl (see also lblimg)
% was computed using the code that accompanies the paper
% "Analyzing spatially-varying blur" by Chakrabarti, Zickler, and Freeman.
% This code, and the image stored in the variable img, are available at
% http://www.eecs.harvard.edu/~ayanc/svblur/

% Here is a Bibtex citation for the paper:
% @inproceedings{chakrabarti2010analyzing,
%   title={Analyzing spatially-varying blur},
%   author={Chakrabarti, A. and Zickler, T. and Freeman, W. T.},
%   booktitle={Computer Vision and Pattern Recognition (CVPR), 2010 IEEE Conference on},
%   pages={2512--2519},
%   year={2010},
%   organization={IEEE}
% }

[numRows,numCols,check] = size(img);

gamma = 1/3000;

U = zeros(numRows,numCols,2);
U(:,:,1) = lbl;
U(:,:,2) = 1 - lbl;


nrMask = 7; 
ncMask = 7;
kernels = zeros(nrMask,ncMask,2);
kernels(4,:,1) = (1/7)*[1 1 1 1 1 1 1];
kernels(4,:,2) = [0 0 0 1 0 0 0];


img = double(img);

mn = min(img(:));
I = img - mn;
mx = max(I(:));
I = I/mx;


% Deblur image using primal-dual Douglas-Rachford.    
t = .003;
beta = 20000;
overRelax = 1.9;

paramsPDDR = struct('maxIter',maxIter,'showTrigger',showTrigger,'t',t,'beta',beta,...
                'gamma',gamma,'overRelax',overRelax);  

tic
[xPDDR1,costsPDDR1] = tvDeblur_varBlur_freeBCs_DR(I(:,:,1),kernels,U,paramsPDDR);  
[xPDDR2,costsPDDR2] = tvDeblur_varBlur_freeBCs_DR(I(:,:,2),kernels,U,paramsPDDR);   
[xPDDR3,costsPDDR3] = tvDeblur_varBlur_freeBCs_DR(I(:,:,3),kernels,U,paramsPDDR);    
timeDeblur = toc;

xPDDR = cat(3,xPDDR1,xPDDR2,xPDDR3);
        
figure('Name','deblurred by Primal-dual Douglas-Rachford')
imshow(xPDDR,[])
