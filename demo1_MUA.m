%% demo1 MUA
% =========================================================================
%
% Demo of example 1 (DC1) of the paper.
%   
% [1] Borsoi, R. A., Imbiriba, T., Bermudez, J. C. M., Richard, C. A Fast Multiscale Spatial Regularization for Sparse Hyperspectral Unmixing. IEEE Geoscience and Remote Sensing Letters, 2018.
%   
%   Sparse unmixing using a multiscale spatial regularization method based on the 
%   segmentation/superpixel decomposition of the HSI. Other techniques
%   are included for comparisson (SUNSAL, SUNSAL-TV and S2WSP). Please refer to [1]
%   for more details.
% 
% 
% Original code by Iordache and Bioucas-Dias
% Modified by: Ricardo Borsoi, 2017
% =========================================================================
%
% 
% 
% This demo illustrates the sunsal_TV sparse regression algorithm
% introduced in the paper 
%
%  M.-D. Iordache, J. Bioucas-Dias, and A. Plaza, "Total variation spatial 
%  regularization for sparse hyperspectral unmixing", IEEE Transactions on 
%  Geoscience and Remote Sensing, vol. PP, no. 99, pp. 1-19, 2012.
%
% which solves the optimization problem
%
%   min  0.5*||AX-Y||^2_F + lambda_1 ||X||_{1,1} + lambda_tv TV(X) 
%   X>=0
%
%
%  Demo parameters:
%     p = 5                             % number of endmembers
%     SNR = 40 dB      
%     size(A) = [220, 240]              % size of the library
%     min angle(a_i, a_j) = 4.44 degs   % minimum angle between any two
%                                       % elements of A
%       
%  Notes:
%
%    You may change the demo parameters, namely SNR, the noise correlation,
%    the size of dictionary A by changing min_angle, and the true endmember 
%    matrix M, which, in any case, must contain p=5 columns. 
% 
%   Please keep in mind the following:
%
%     a) sunsal  adapts automatically  the ADMM parameter for 
%        convergence speed 
%  
%     b) sunsal_tv deoes not adapts automatically  the ADMM parameter. 
%        So the inputted parameter mu has a  critical impact on the
%        convergence speed
%
%     c) the regularization parameters  were hand tuned for optimal
%        performance.
% 
% Author: Jose Bioucas Dias, August 2012
%


close all
clear all
clc

mkdir('examples/DC0')
% mkdir('mat_data')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of end members
p = 5;  % fixed for this demo

%SNR in dB
SNR = 30; 20;30;30;20;
% noise bandwidth in pixels of the noise  low pass filter (Gaussian)
bandwidth = 10000; % 10000 == iid noise
%bandwidth = 5*pi/224; % colored noise 


% define random states
rand('state',10);
randn('state',10);


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gererate fractional abundances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pure pixels
x1 = eye(p);

% mixtures with two materials
x2 = x1 + circshift(eye(p),[1 0]);

% mixtures with three materials
x3 = x2 + circshift(eye(p),[2 0]);

% mixtures with four  materials
x4 = x3 + circshift(eye(p),[3 0]);

% mixtures with four  materials
x5 = x4 + circshift(eye(p),[4 0]);


% normalize
x2 = x2/2;
x3 = x3/3;
x4 = x4/4;
x5 = x5/5;


% background (random mixture)
%x6 = dirichlet(ones(p,1),1)';
x6 = [0.1149 0.0741  0.2003 0.2055, 0.4051]';   % as in the paper

% build a matrix
xt = [x1 x2 x3 x4 x5 x6];


% build image of indices to xt
imp = zeros(3);
imp(2,2)=1;


imind = [imp*1  imp*2 imp* 3 imp*4 imp*5;
    imp*6  imp*7 imp* 8 imp*9 imp*10;
    imp*11  imp*12 imp*13 imp*14 imp*15;
    imp*16  imp*17 imp* 18 imp*19 imp*20;
    imp*21  imp*22 imp* 23 imp*24 imp*25];

imind = kron(imind,ones(5));

% set backround index
imind(imind == 0) = 26;

% generare frectional abundances for all pixels
[nl,nc] = size(imind);
np = nl*nc;     % number of pixels
for i=1:np
    X(:,i) = xt(:,imind(i));
end


% reshape and plot some figures
Xim = reshape(X',nl,nc,p);

%  image endmember 1
figure(1)
imagesc(Xim(:,:,5))
title('Frational abundance of endmember 5')

%  Size of the image
nl = size(imind,1);
nc = size(imind,2);


%% Plot all endmembers
figure
subplot(1,5,1)
imagesc(Xim(:,:,1), [0 1])
axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])

subplot(1,5,2)
imagesc(Xim(:,:,2), [0 1])
axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])

subplot(1,5,3)
imagesc(Xim(:,:,3), [0 1])
axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])

subplot(1,5,4)
imagesc(Xim(:,:,4), [0 1])
axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])

h2 = subplot(1,5,5);
imagesc(Xim(:,:,5), [0 1])
axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
originalSize2 = get(gca, 'Position');

h=colorbar; 
set(h2, 'Position', originalSize2);
set(h,'fontsize',5);

colormap jet

print('examples/DC0/true_abundances_DC0','-depsc')
close all

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% buid the dictionary 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load USGS_1995_Library.mat
%  order bands by increasing wavelength
[dummy index] = sort(datalib(:,1));
A =  datalib(index,4:end);
names = names(4:end,:);


% prune the library 
% min angle (in degres) between any two signatures 
% the larger min_angle the easier is the sparse regression problem
min_angle = 4.44;       
[A, index] = prune_library2(A,min_angle); % 240  signature 
names = names(index',:);

% order  the columns of A by decreasing angles 
[A, index, angles] = sort_library_by_angle(A);
names = names(index',:);


% Names of the first 10 ordered materials, with 4.44 deg. prunning:
% 1 - Jarosite GDS99 K,Sy 200C
% 2 - Jarosite GDS101 Na,Sy 200
% 3 - Anorthite HS349.3B 
% 4 - Calcite WS272 
% 5 - Alunite GDS83 Na63 
% 6 - Howlite GDS155
% 7 - Corrensite CorWa-1
% 8 - Fassaite HS118.3B  
% 9 - Adularia GDS57 Orthoclase  
% 10 - Andradite NMNH113829 


%% select p endmembers  from A

% angles (a_1,a_j) \sisizemeq min_angle)
% supp = 1:p;
supp = [2 3 4 5 6]; % dont get two Jarosites

% % Sample endmembers at random
% supp = randsample(size(A,2), p);

M = A(:,supp);
[L,p] = size(M);  % L = number of bands; p = number of material


%%
%---------------------------------
% generate  the observed  data X
%---------------------------------

% set noise standard deviation
sigma = sqrt(sum(sum((M*X).^2))/np/L/10^(SNR/10));
% generate Gaussian iid noise
noise = sigma*randn(L,np);


% make noise correlated by low pass filtering
% low pass filter (Gaussian)
filter_coef = exp(-(0:L-1).^2/2/bandwidth.^2)';
scale = sqrt(L/sum(filter_coef.^2));
filter_coef = scale*filter_coef;
noise = idct(dct(noise).*repmat(filter_coef,1,np));

%  observed spectral vector
Y = M*X + noise;


% create  true X wrt  the library A
n = size(A,2);
N = nl*nc;
XT = zeros(n,N);
XT(supp,:) = X;







%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularization using K-means
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define random states
rand('state',10);
randn('state',10);

% Select parameters acording to SNR
if SNR == 40
    % 40db
    error('No parameters for this SNR!')
elseif SNR == 30
    lambda1_sp = 0.005;
    lambda2_sp = 0.05;
    beta       = 10; % opt = 30, flat surface
    kmensCsize = 7;
    
elseif SNR == 20
    lambda1_sp = 0.005;
    lambda2_sp = 0.5;
    beta       = 30;
    kmensCsize = 13;
end


    
rmpath('vlfeat-0.9.20/toolbox/noprefix')

% Determine the number of clusters
NclustKmeans = floor(N/kmensCsize^2);



tic

Y2 = reshape(Y', nl, nc, L);   
Y2a = Y2;
% reorder and rescale data into 2-D array
[numRows,numCols,numSpectra] = size(Y2);
scfact = mean(reshape(sqrt(sum(Y2.^2,3)), numRows*numCols, 1));
Y2 = Y2./scfact;
imgVec = reshape(Y2, [numRows*numCols numSpectra]);

% apply k-means
[IDX, C] = kmeans(Y', NclustKmeans, 'Distance', 'correlation','Start',rand(NclustKmeans,L));

IDX = IDX - 1;
spSegs = reshape(IDX', nl, nc, 1); 
numSuperpixels = NclustKmeans;


% ------
% Unmix the clusters

Y3 = zeros(size(Y2));
avg_superpx = zeros(1, numSuperpixels+1, L);

for i=0:numSuperpixels
    [rowi, coli] = find(spSegs==i);
    
    for j=1:length(rowi)
        % Averages all pixels inside each superpixel
        if j == 1
            avg_superpx(1,i+1,:) = (1/length(rowi)) * Y2a(rowi(j),coli(j),:);
        else
            avg_superpx(1,i+1,:) = avg_superpx(1,i+1,:) + (1/length(rowi)) * Y2a(rowi(j),coli(j),:);
        end
    end
    
    % This is optional (for visualization)
    for j=1:length(rowi)
        Y3(rowi(j),coli(j),:) = avg_superpx(1,i+1,:);
    end
end



% %%
% Unmix each superpixel individually
[X_hat_l1_t_kmeans] = sunsal(A,squeeze(avg_superpx)','lambda',lambda1_sp,'ADDONE','no','POSITIVITY','yes', ...
                       'TOL',1e-4, 'AL_iters',2000,'verbose','yes');

% Re-attribute the abundances for the entire matrix
temp = zeros(size(Y2,1), size(Y2,2), n);
for i=0:numSuperpixels
    [rowi, coli] = find(spSegs==i);
    
    % Attributes unmixing result to all pixels in a voxel
    for j=1:length(rowi)
        temp(rowi(j),coli(j),:) = X_hat_l1_t_kmeans(:,i+1);
    end
end

X_hat_l1_kmeans = reshape(temp, [size(Y2,1)*size(Y2,2) n])';

% constrained least squares l2-l1 
[X_hat_l1_kmeans] =  sunsal_spreg(A,Y,X_hat_l1_kmeans,beta,'lambda',lambda2_sp,'ADDONE','no','POSITIVITY','yes', ...
                    'TOL',1e-4, 'AL_iters',2000,'verbose','yes');
timeKmeans = toc;

SRE_l1_kmeans = 20*log10(norm(XT,'fro')/norm(X_hat_l1_kmeans-XT,'fro'));
disp(SRE_l1_kmeans)
















%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiscale regularization using segmentation (watershed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('HSI_segmentation'))


% Select parameters acording to SNR, considering dirichlet data
if SNR == 40
    % 40db
    error('No parameters for this SNR!')
elseif SNR == 30
    % 30db
    lambda1_sp = 0.005;
    lambda2_sp = 0.1;
    beta       = 30; % opt = 30, flat surface
    sideBPT    = 10;
elseif SNR == 20
    % 20db
    lambda1_sp = 0.005;
    lambda2_sp = 0.1;
    beta       = 30;
    sideBPT    = 14;
end


% Set thenumber of partitions
numBPTsegs = floor(N/sideBPT^2);



tic

Y2 = reshape(Y', nl, nc, L);   
Y2a = Y2;
% reorder and rescale data into 2-D array
[numRows,numCols,numSpectra] = size(Y2);
scfact = mean(reshape(sqrt(sum(Y2.^2,3)), numRows*numCols, 1));
Y2 = Y2./scfact;
imgVec = reshape(Y2, [numRows*numCols numSpectra]);


% Main for BPT construction and processing
% Initial segmentation

% It is assumed that a N-D image "im" has been already loaded, and a RGB
% version "imrgb" of this image is also available (for segmentation display
% purposes)

% first get a pre-segmentation of the image (optional, but advised...)
% Several methods are possible, such as watershed segmentation, SLIC or
% mean shift clustering (code is available online for SLIC and MSC)
% Here is provided the multidimensional watershed algorithm
% Best pre-segmentation algorithm is the mean shift clustering without any
% doubt

im = Y2a;
imrgb = Y2a(:,:,[10 20 90]);

% Segmentation by multidimensional watershed
% ------------------------------------------
SEG = multidimwatershed(im,'supremum');
initsegmap = SEG.whed;

% build data structures

% The two available region models are the region-wise mean vector "R_mean" 
% and the region-wise collection of histograms "R_hist", with a bunch of
% associated metrics (Euclidean, L1 and Linfinity norms, SAM,
% Kullback-Leibleir and Jensen-Shannon for the mean model, Battacharyya and
% Diffusion distance  for the histogram model)

% The field specmerging in structure DAT has to be specified as
% "merging_mean" or "merging_hist" according to the chosen region model
DAT = struct;
DAT.data = im;
DAT.initsegmap = initsegmap;
DAT.regionmodel = @R_mean;
DAT.mergingcriterion = @O_SAM;
DAT.specmerging = @merging_mean;
DAT.prioritysize = @priority_size;

% build BPT
T = buildBPT(DAT);

% Process/prune BPT
% BPT processing can be done in different way, and must be specifically
% tuned to the desired goal.

% A simple way to process the BPT is to extract a given number of regions N,
% which will constitute a partition featuring the N most dissimilar regions
% created during the construction of the BPT.

prunedtree = pruneBPTnbregions(T,numBPTsegs); % N=10


% Display partition
segregions = retrievesegmentation(prunedtree,initsegmap,'regions',imrgb);
% drawborders(imrgb,segregions,'red');
% imrgbmean = displaysegmentationfalsecolors(segregions,imrgb);
% figure, imshow(imrgbmean)


temppp = -ones(size(segregions));
tempVals = sort(unique(segregions));
for iii = 1:length(tempVals)
    temppp(segregions == tempVals(iii)) = iii;
end
temppp = temppp-1;

% figure, imagesc(segregions)
% figure, imagesc(temppp)

if any(temppp == -1)
    error('Error in the segmentation regions partitioning!')
end



spSegs = temppp;
numSuperpixels = numBPTsegs;

% % compute superpixels
% disp('Computing SLIC Superpixels...');
% spSegs = vl_slic(single(Y2), slic_size, slic_reg);
% numSuperpixels = double(max(spSegs(:)))+1; 

% ------
% Unmix the clusters

Y3 = zeros(size(Y2));
avg_superpx = zeros(1, numSuperpixels+1, L);

for i=0:numSuperpixels
    [rowi, coli] = find(spSegs==i);
    
    for j=1:length(rowi)
        % Averages all pixels inside each superpixel
        if j == 1
            avg_superpx(1,i+1,:) = (1/length(rowi)) * Y2a(rowi(j),coli(j),:);
        else
            avg_superpx(1,i+1,:) = avg_superpx(1,i+1,:) + (1/length(rowi)) * Y2a(rowi(j),coli(j),:);
        end
    end
    
    % This is optional (for visualization)
    for j=1:length(rowi)
        Y3(rowi(j),coli(j),:) = avg_superpx(1,i+1,:);
    end
end


% %%
% Unmix each superpixel individually
[X_hat_l1_t_BPT] = sunsal(A,squeeze(avg_superpx)','lambda',lambda1_sp,'ADDONE','no','POSITIVITY','yes', ...
                       'TOL',1e-4, 'AL_iters',2000,'verbose','yes');

% Re-attribute the abundances for the entire matrix
temp = zeros(size(Y2,1), size(Y2,2), n);
for i=0:numSuperpixels
    [rowi, coli] = find(spSegs==i);
    
    % Attributes unmixing result to all pixels in a voxel
    for j=1:length(rowi)
        temp(rowi(j),coli(j),:) = X_hat_l1_t_BPT(:,i+1);
    end
end

X_hat_l1_BPT = reshape(temp, [size(Y2,1)*size(Y2,2) n])';

% constrained least squares l2-l1 
[X_hat_l1_BPT] =  sunsal_spreg(A,Y,X_hat_l1_BPT,beta,'lambda',lambda2_sp,'ADDONE','no','POSITIVITY','yes', ...
                    'TOL',1e-4, 'AL_iters',2000,'verbose','yes');
timeBPT = toc;

SRE_l1_BPT = 20*log10(norm(XT,'fro')/norm(X_hat_l1_BPT-XT,'fro'));
disp(SRE_l1_BPT)

% % endmember no. 5
% X_hat_l1_spreg_im = reshape(X_hat_l1_spreg', nl,nc,n);       
% figure, imagesc(X_hat_l1_spreg_im(:,:,supp(5)))
% title('Spreg - Frational abundance of endmember 5')








%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiscale regularization using superpixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Decompose the image into homogeneous regions and apply clustering
% individually 

addpath vlfeat-0.9.20
run('vlfeat-0.9.20/toolbox/vl_setup')


% Select parameters acording to SNR, considering dirichlet data
if SNR == 40
    % 40db
    lambda1_sp = 1e-3;
    lambda2_sp = 0.01;
    beta       = 10; % opt = 30, flat surface
    slic_size  = 5;
    slic_reg   = 0.01;
    
elseif SNR == 30
    % 30db
    lambda1_sp = 0.007;
    lambda2_sp = 0.05;
    beta       = 10; % opt = 30, flat surface
    slic_size  = 5;
    slic_reg   = 0.005;
    
elseif SNR == 20
    % 20db
    lambda1_sp = 0.03;
    lambda2_sp = 0.1;
    beta       = 30;
    slic_size  = 6;
    slic_reg   = 0.005;
end


Y2 = reshape(Y', nl, nc, L);   
Y2a = Y2;

tic

% reorder and rescale data into 2-D array
[numRows,numCols,numSpectra] = size(Y2);
scfact = mean(reshape(sqrt(sum(Y2.^2,3)), numRows*numCols, 1));
Y2 = Y2./scfact;
imgVec = reshape(Y2, [numRows*numCols numSpectra]);

% compute superpixels
disp('Computing SLIC Superpixels...');
spSegs = vl_slic(single(Y2), slic_size, slic_reg);
numSuperpixels = double(max(spSegs(:)))+1; 


% ------
% Unmix the superpixels

Y3 = zeros(size(Y2));
avg_superpx = zeros(1, numSuperpixels+1, L);

for i=0:numSuperpixels
    [rowi, coli] = find(spSegs==i);
    
    for j=1:length(rowi)
        % Averages all pixels inside each superpixel
        if j == 1
            avg_superpx(1,i+1,:) = (1/length(rowi)) * Y2a(rowi(j),coli(j),:);
        else
            avg_superpx(1,i+1,:) = avg_superpx(1,i+1,:) + (1/length(rowi)) * Y2a(rowi(j),coli(j),:);
        end
    end
    
    % This is optional (for visualization)
    for j=1:length(rowi)
        Y3(rowi(j),coli(j),:) = avg_superpx(1,i+1,:);
    end
end


% %%
% Unmix each superpixel individually
[X_hat_l1_suppx] = sunsal(A,squeeze(avg_superpx)','lambda',lambda1_sp,'ADDONE','no','POSITIVITY','yes', ...
                       'TOL',1e-4, 'AL_iters',2000,'verbose','yes');

% Re-attribute the abundances for the entire matrix
temp = zeros(size(Y2,1), size(Y2,2), n);
for i=0:numSuperpixels
    [rowi, coli] = find(spSegs==i);
    
    % Attributes unmixing result to all pixels in a voxel
    for j=1:length(rowi)
        temp(rowi(j),coli(j),:) = X_hat_l1_suppx(:,i+1);
    end
end

X_hat_l1_spreg = reshape(temp, [size(Y2,1)*size(Y2,2) n])';

% constrained least squares l2-l1 
[X_hat_l1_spreg] =  sunsal_spreg(A,Y,X_hat_l1_spreg,beta,'lambda',lambda2_sp,'ADDONE','no','POSITIVITY','yes', ...
                    'TOL',1e-4, 'AL_iters',2000,'verbose','yes');
timeMscale = toc;

SRE_l1_spreg = 20*log10(norm(XT,'fro')/norm(X_hat_l1_spreg-XT,'fro'));


% % endmember no. 5
% X_hat_l1_spreg_im = reshape(X_hat_l1_spreg', nl,nc,n);       
% figure, imagesc(X_hat_l1_spreg_im(:,:,supp(5)))
% title('Spreg - Frational abundance of endmember 5')






% ------------------------------------------
% display image of superpixels (optional)
[sx,sy] = vl_grad(double(spSegs), 'type', 'forward') ;
s = find(sx | sy) ;
imgColor = Y2(:,:,[29 15 12]);
imgColor = uint8(255*(imgColor - min(imgColor(:)))./(max(imgColor(:))-min(imgColor(:))));
imgS = imgColor; 
imgS([s s+numel(imgColor(:,:,1)) s+2*numel(imgColor(:,:,1))]) = 0;
figure; imshow(imgS);



















%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUNSAL and SUNSAL_TV solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                
%% constrained least squares l2-l1   (SUnSAL)              
% Select parameters acording to SNR, considering dirichlet data
if SNR == 40
    % 40db
    lambda = 0.01;
    
elseif SNR == 30
    % 30db
    lambda = 0.1;
    
elseif SNR == 20
    % 20db
    lambda = 0.7;
end

tic
[X_hat_l1] =  sunsal(A,Y,'lambda',lambda,'ADDONE','no','POSITIVITY','yes', ...
                    'TOL',1e-4, 'AL_iters',2000,'verbose','yes');
timeSunsal = toc;

SRE_l1 = 20*log10(norm(XT,'fro')/norm(X_hat_l1-XT,'fro'));



%% constrained least squares l2-l1-TV   (SUnSAL-TV)      
% Select parameters acording to SNR, considering dirichlet data
if SNR == 40
    % 40db
    lambda    = 0.001;
    lambda_TV = 0.003;
    
elseif SNR == 30
    % 30db
    lambda    = 0.007;
    lambda_TV = 0.01;
    
elseif SNR == 20
    % 20db
    lambda    = 0.05;
    lambda_TV = 0.05;
end



tic
[X_hat_tv,res,rmse] = sunsal_tv(A,Y,'MU',0.05,'POSITIVITY','yes','ADDONE','no', ...
                              'LAMBDA_1',lambda,'LAMBDA_TV', lambda_TV, 'TV_TYPE','niso',...
                              'IM_SIZE',[nl,nc],'AL_ITERS',200, 'TRUE_X', XT,  'VERBOSE','yes');
timeTV = toc;

SRE_tv = 20*log10(norm(XT,'fro')/norm(X_hat_tv-XT,'fro')); 
             




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constrained least squares using l2-l1-swSp solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
% Select parameters acording to SNR
if SNR == 40
    % 40db
    
elseif SNR == 30
    % 30db
    lambda_swsp = 5e-3;
    
elseif SNR == 20
    % 20db
    lambda_swsp = 100e-3;
end


tic
[X_hat_l21LC,res2,rmse_ni] = sunsal_tv_lw_sp(A,Y,'MU',0.5,'POSITIVITY','yes','ADDONE','no', ...
                               'LAMBDA_1',lambda_swsp,'IM_SIZE',[nl,nc],'AL_ITERS',5,'TRUE_X', XT, 'VERBOSE','yes');  
time_swsp = toc;

SRE_swsp = 20*log10(norm(XT,'fro')/norm(X_hat_l21LC-XT,'fro')); 
   


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('\n\n SIGNAL-TO-RECONSTRUCTION ERRORS (SRE)\n\n')

fprintf('SRE-l1 = %2.3f\n', SRE_l1)
fprintf('SRE-TV = %2.3f\n', SRE_tv)
fprintf('SRE_SPx-Mscale = %2.3f\n', SRE_l1_spreg)
fprintf('SRE_kmeans = %2.3f\n', SRE_l1_kmeans)
fprintf('SRE_BPT = %2.3f\n', SRE_l1_BPT)
fprintf('SRE_swsp = %2.3f\n\n', SRE_swsp)


fprintf('\n\n EXECUTION TIME \n\n')

fprintf('TIME-l1 = %2.3f\n', timeSunsal)
fprintf('TIME-TV = %2.3f\n', timeTV)
fprintf('TIME_SPx-Mscale = %2.3f\n', timeMscale)
fprintf('TIME-kmeans = %2.3f\n', timeKmeans)
fprintf('TIME-BPT = %2.3f\n', timeBPT)
fprintf('TIME-SWSP = %2.3f\n\n', time_swsp)  

        





%%


X_hat_l1_im       = reshape(X_hat_l1', nl,nc,n);  
X_hat_tv_im       = reshape(X_hat_tv', nl,nc,n);  
X_hat_l1_spreg_im = reshape(X_hat_l1_spreg', nl,nc,n);  
X_hat_kmeans_im   = reshape(X_hat_l1_kmeans', nl,nc,n);  
X_hat_segment_im  = reshape(X_hat_l1_BPT', nl,nc,n);  
X_hat_l21LC_im    = reshape(X_hat_l21LC', nl,nc,n); 


em_idx = 2;

figure
% [ha, pos] = tight_subplot(1,5,[.025 .025],[.05 .05],[.05 .15]);
[ha, pos] = tight_subplot(1,5,[.025 .015],[.05 .05],[.05 .15]);


axes(ha(1));
imagesc(Xim(:,:,(em_idx)), [0 1])
axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])

axes(ha(2));
imagesc(X_hat_tv_im(:,:,supp(em_idx)), [0 1])
axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])

axes(ha(3));
imagesc(X_hat_l21LC_im(:,:,supp(em_idx)), [0 1])
axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])

axes(ha(4));
imagesc(X_hat_segment_im(:,:,supp(em_idx)), [0 1])
axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])

axes(ha(5));
imagesc(X_hat_l1_spreg_im(:,:,supp(em_idx)), [0 1])
axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
originalSize2 = get(gca, 'Position');

h=colorbar; 
set(ha(5), 'Position', originalSize2);
% set(h,'fontsize',5);
set(h,'fontsize',8);

colormap jet
print(strcat('examples/DC1b/estim_abundances_DC0_SNR',num2str(SNR),'_tght'),'-dpdf')


























