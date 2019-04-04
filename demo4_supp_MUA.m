%% supplemental material, demo4 MUA
% =========================================================================
%
% Demo of example 4 (Samson image) of the paper:
% 
% [1] Borsoi, R. A., Imbiriba, T., Bermudez, J. C. M., Richard, C. A Fast Multiscale Spatial Regularization for Sparse Hyperspectral Unmixing. IEEE Geoscience and Remote Sensing Letters, 2018.
% 
% This example is discussed in the supplemental material of the extended 
% manuscript, available at: https://arxiv.org/pdf/1712.01770.pdf
%   
%   Sparse unmixing using a multiscale spatial regularization method based on the 
%   segmentation/superpixel decomposition of the HSI. Other techniques
%   are included for comparisson (SUNSAL, SUNSAL-TV and S2WSP). Please refer to [1]
%   for more details.
% 
% 
% Written by: Ricardo Borsoi, 2018
% =========================================================================

close all
clear all
clc

mkdir('examples')
%mkdir('mat_data')

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load hyperspectral image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% HERE YOU NEED TO LOAD THE SAMSON IMAGE FILE.
addpath('real_data')
load('samson_1.mat')


% define random states
rand('state',10);
randn('state',10);


Y = V;
Yim = reshape(Y',nRow,nCol,nBand);

% set constants
ptrue = 3;
nl = nRow;
nc = nCol;
np = nl*nc;
L  = size(Y,1);
N  = nl*nc;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the dictionary 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('spectral_library_samson.mat')
n = size(A,2);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularization using K-means
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define random states
rand('state',10);
randn('state',10);


lambda1_sp = 0.005;
lambda2_sp = 0.01;
beta       = 1; % opt = 30, flat surface
kmensCsize = 8;


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

% run kmeans
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








%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiscale regularization using segmentation (watershed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('HSI_segmentation'))


% parameters -------
lambda1_sp = 0.001;
lambda2_sp = 0.05;
beta       = 1;
sideBPT    = 11;

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





% ----------------------------------------------------------------
% display image of superpixels (optional)
[sx,sy] = vl_grad(double(spSegs), 'type', 'forward') ;
s = find(sx | sy);

imgColor = Y2(:,:,[100 80 50]);
imgColor = uint8(255*(imgColor - min(imgColor(:)))./(max(imgColor(:))-min(imgColor(:))));
imgS = imgColor; 
imgS([s s+numel(imgColor(:,:,1)) s+2*numel(imgColor(:,:,1))]) = 0;

figure, imshow(imgColor)
figure; imshow(1.5*imgS)

















%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiscale regularization using superpixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Decompose the image into homogeneous regions and apply clustering
% individually 

addpath SLIC_DimensionReduction
addpath SSSE
addpath vlfeat-0.9.20
run('vlfeat-0.9.20/toolbox/vl_setup')



% Parameters -------
lambda1_sp = 0.003;
lambda2_sp = 0.03;
beta       = 3;
slic_size  = 7;
slic_reg   = 0.00125;




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



% ---------------------------------------------
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


lambda = 0.01;

tic
[X_hat_l1] =  sunsal(A,Y,'lambda',lambda,'ADDONE','no','POSITIVITY','yes', ...
                    'TOL',1e-4, 'AL_iters',2000,'verbose','yes');
timeSunsal = toc;


                
%% SUnSAL-TV

lambda    = 0.005;
lambda_TV = 0.007;

XT = zeros(n,N);
tic
[X_hat_tv,res,rmse_ni] = sunsal_tv(A,Y,'MU',0.05,'POSITIVITY','yes','ADDONE','no', ...
                               'LAMBDA_1',lambda,'LAMBDA_TV', lambda_TV, 'TV_TYPE','niso',...
                               'IM_SIZE',[nl,nc],'AL_ITERS',200, 'TRUE_X', XT,  'VERBOSE','yes');
timeTV = toc;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constrained least squares using l2-l1-swSp solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda_swsp = 10e-3;

XT = zeros(n,N);
tic
[X_hat_l21LC,res2,rmse_ni] = sunsal_tv_lw_sp(A,Y,'MU',0.5,'POSITIVITY','yes','ADDONE','no', ...
                               'LAMBDA_1',lambda_swsp,'IM_SIZE',[nl,nc],'AL_ITERS',5,'TRUE_X', XT, 'VERBOSE','yes');  
time_swsp = toc;






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('\n\n EXECUTION TIME \n\n')

fprintf('TIME-l1 = %2.3f\n', timeSunsal)
fprintf('TIME-TV = %2.3f\n', timeTV)
fprintf('TIME_SPx-Mscale = %2.3f\n', timeMscale)
fprintf('TIME-kmeans = %2.3f\n', timeKmeans)
fprintf('TIME-BPT = %2.3f\n', timeBPT)
fprintf('TIME-SWSP = %2.3f\n\n', time_swsp)  

        
        



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot abundance maps for the desired materials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nLib1 = size(lib1,2);
nLib2 = size(lib2,2);
nLib3 = size(lib3,2);
ptrue = 3;

X_hat_l12 = [mean(X_hat_l1(1:nLib1,:),1);...
             mean(X_hat_l1((nLib1+1):(nLib1+nLib2),:),1); ...
             mean(X_hat_l1((nLib1+nLib2+1):(nLib1+nLib2+nLib3),:),1)];
         
X_hat_tv2 = [mean(X_hat_tv(1:nLib1,:),1);...
             mean(X_hat_tv((nLib1+1):(nLib1+nLib2),:),1); ...
             mean(X_hat_tv((nLib1+nLib2+1):(nLib1+nLib2+nLib3),:),1)];
         
X_hat_l1_spreg2 = [mean(X_hat_l1_spreg(1:nLib1,:),1);...
                   mean(X_hat_l1_spreg((nLib1+1):(nLib1+nLib2),:),1); ...
                   mean(X_hat_l1_spreg((nLib1+nLib2+1):(nLib1+nLib2+nLib3),:),1)];
         
X_hat_l1_kmeans2 = [mean(X_hat_l1_kmeans(1:nLib1,:),1);...
                    mean(X_hat_l1_kmeans((nLib1+1):(nLib1+nLib2),:),1); ...
                    mean(X_hat_l1_kmeans((nLib1+nLib2+1):(nLib1+nLib2+nLib3),:),1)];
         
X_hat_l1_BPT2 = [mean(X_hat_l1_BPT(1:nLib1,:),1);...
                 mean(X_hat_l1_BPT((nLib1+1):(nLib1+nLib2),:),1); ...
                 mean(X_hat_l1_BPT((nLib1+nLib2+1):(nLib1+nLib2+nLib3),:),1)];

X_hat_l21LC2 = [mean(X_hat_l21LC(1:nLib1,:),1);...
                mean(X_hat_l21LC((nLib1+1):(nLib1+nLib2),:),1); ...
                mean(X_hat_l21LC((nLib1+nLib2+1):(nLib1+nLib2+nLib3),:),1)];
        
             
             
             
% also need to normalize
X_hat_l12 = X_hat_l12 ./repmat(sum(X_hat_l12), [ptrue 1]);
X_hat_tv2 = X_hat_tv2 ./repmat(sum(X_hat_tv2), [ptrue 1]);
X_hat_l1_spreg2 = X_hat_l1_spreg2 ./repmat(sum(X_hat_l1_spreg2), [ptrue 1]);
X_hat_l1_kmeans2 = X_hat_l1_kmeans2 ./repmat(sum(X_hat_l1_kmeans2), [ptrue 1]);
X_hat_l1_BPT2 = X_hat_l1_BPT2 ./repmat(sum(X_hat_l1_BPT2), [ptrue 1]);
X_hat_l21LC2 = X_hat_l21LC2 ./repmat(sum(X_hat_l21LC2), [ptrue 1]);


% reshape as image
X_hat_l1_im2 = reshape(X_hat_l12', nl,nc,ptrue);  
X_hat_tv_im2 = reshape(X_hat_tv2', nl,nc,ptrue);  
X_hat_l1_spreg_im2 = reshape(X_hat_l1_spreg2', nl,nc,ptrue);   
X_hat_kmeans_im2 = reshape(X_hat_l1_kmeans2', nl,nc,ptrue);  
X_hat_segment_im2 = reshape(X_hat_l1_BPT2', nl,nc,ptrue);  
X_hat_l21LC2_im2 = reshape(X_hat_l21LC2', nl,nc,ptrue);  




figure
[ha, pos] = tight_subplot(3,4,[.025 .025],[.05 .05],[.05 .15]);

numAlgs = 4; % numebr of algorithms

for i=1:ptrue
%     axes(ha((i-1)*numAlgs + 1));
%     imagesc(Xim(:,:,i), [0 1]), axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
    
%     axes(ha((i-1)*numAlgs + 2));
%     imagesc(X_hat_tv_im2(:,:,i), [0 1]), axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])

    axes(ha((i-1)*numAlgs + 1));
    imagesc(X_hat_tv_im2(:,:,i), [0 1]), axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
    
    axes(ha((i-1)*numAlgs + 2));
    imagesc(X_hat_l21LC2_im2(:,:,i), [0 1]), axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])

    axes(ha((i-1)*numAlgs + 3));
    imagesc(X_hat_segment_im2(:,:,i), [0 1]), axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
    
    axes(ha((i-1)*numAlgs + 4));
    imagesc(X_hat_l1_spreg_im2(:,:,i), [0 1]), axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
    originalSize2 = get(gca, 'Position');
    
    h=colorbar; 
    set(ha((i-1)*numAlgs + 4), 'Position', originalSize2);
%     set(h,'fontsize',5);
    set(h,'fontsize',8);
end


axes(ha((1-1)*numAlgs + 1)); ylabel('Soil','interpreter','latex')
axes(ha((2-1)*numAlgs + 1)); ylabel('Tree','interpreter','latex')
axes(ha((3-1)*numAlgs + 1)); ylabel('Water','interpreter','latex')


axes(ha(end-3)); xlabel('SUnSAL-TV','interpreter','latex')
axes(ha(end-2)); xlabel('S$^2$WSU','interpreter','latex')
axes(ha(end-1)); xlabel('MUA (BPT)','interpreter','latex')
axes(ha(end)); xlabel('MUA (SLIC)','interpreter','latex')




%%

mkdir('examples')
print('examples/estim_abundances_samson','-dpdf')


figure, imagesc(5*Yim(:,:,[70 50 35]))
set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
print('examples/samson_im','-dpdf')









