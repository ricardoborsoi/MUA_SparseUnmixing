%% demo3 MUA
% =========================================================================
%
% Demo of example 3 (Cuprite image) of the paper:
% 
% [1] Borsoi, R. A., Imbiriba, T., Bermudez, J. C. M., Richard, C. A Fast Multiscale Spatial Regularization for Sparse Hyperspectral Unmixing. IEEE Geoscience and Remote Sensing Letters, 2018.
% 
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

mkdir('examples/realImg')
% mkdir('mat_data')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load hyperspectral image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% HERE YOU NEED TO LOAD THE CUPRITE IMAGE FILE. YOU CAN DOWNLOAD IT FROM http://www.lx.it.pt/~bioucas/code.htm
addpath real_data
load cuprite_ref.mat


% define random states
rand('state',10);
randn('state',10);

% L - number of bands
nl = Lines;
nc = Columns;
Yim = reshape(x',nl,nc,L);

imagesc(Yim(:,:,50))

nl = size(Yim,1);
nc = size(Yim,2);


% reorder in the form of a 2D matrix
Y = reshape(Yim, [size(Yim,1)*size(Yim,2) L])';


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% buid the dictionary 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load USGS_1995_Library.mat
%  order bands by increasing wavelength
[dummy index] = sort(datalib(:,1));
A =  datalib(index,4:end);
names = names(4:end,:);

% % prune the library 
% % min angle (in degres) between any two signatures 
% % the larger min_angle the easier is the sparse regression problem
% min_angle = 4.44;       
% [A, index] = prune_library2(A,min_angle); % 240  signature 
% names = names(index',:);

% order  the columns of A by decreasing angles 
[A, index, angles] = sort_library_by_angle(A);
names = names(index',:);
namesStr = char(names);

A = A(BANDS,:);

n = size(A,2);




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularization using K-means
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 250*191;

% define random states
rand('state',10);
randn('state',10);


lambda1_sp = 0.001
lambda2_sp = 0.001
beta       = 3
kmensCsize = 11;


    
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

% run k-means
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




% endmember no. 5
X_hat_l1_kmeans_im = reshape(X_hat_l1_kmeans', nl,nc,n);       
figure, imagesc(X_hat_l1_kmeans_im(:,:,420),[0 0.38])
colormap jet
title('Kmeans - Frational abundance of endmember 5')



  










%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiscale regularization using segmentation (watershed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('HSI_segmentation'))


N = 250*191;

% Parameters ------------
lambda1_sp = 0.001;
lambda2_sp = 0.001;
beta       = 3;
sideBPT    = 7;







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



% endmember no. 5
X_hat_l1_BPT_im = reshape(X_hat_l1_BPT', nl,nc,n);       
figure, imagesc(X_hat_l1_BPT_im(:,:,420),[0 0.38])
colormap jet
title('BPT - Frational abundance of endmember 5')





% ----------------------------------------------------------------
% display image of superpixels (optional)
[sx,sy] = vl_grad(double(spSegs), 'type', 'forward') ;
s = find(sx | sy);

imgColor = Y2(:,:,[100 80 50]);
imgColor = uint8(255*(imgColor - min(imgColor(:)))./(max(imgColor(:))-min(imgColor(:))));
imgS = imgColor; 
imgS([s s+numel(imgColor(:,:,1)) s+2*numel(imgColor(:,:,1))]) = 0;

figure, imshow(imgColor)
print('examples/realImg/img_example2','-depsc')
figure; imshow(1.5*imgS)
print('examples/realImg/img_example_BPT','-depsc')
close all










%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiscale regularization using superpixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath SLIC_DimensionReduction
addpath SSSE
addpath vlfeat-0.9.20
run('vlfeat-0.9.20/toolbox/vl_setup')




% Parameters ------------
lambda1_sp = 0.001;
lambda2_sp = 0.001;
beta       = 3;
slic_size  = 5;
slic_reg   = 0.0025;




% Reorder images
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


% ---------------
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



% ----------------------------------------------------------------
% display image of superpixels (optional)
[sx,sy] = vl_grad(double(spSegs), 'type', 'forward') ;
s = find(sx | sy);

imgColor = Y2(:,:,[100 80 50]);
imgColor = uint8(255*(imgColor - min(imgColor(:)))./(max(imgColor(:))-min(imgColor(:))));
imgS = imgColor; 
imgS([s s+numel(imgColor(:,:,1)) s+2*numel(imgColor(:,:,1))]) = 0;

figure, imshow(imgColor)
print('examples/realImg/img_example','-depsc')
figure; imshow(1.5*imgS)
print('examples/realImg/img_example_sppx','-depsc')
close all





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUNSAL and SUNSAL_TV solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = 1e-3;

tic
[X_hat_l1] =  sunsal(A,Y,'lambda',lambda,'ADDONE','no','POSITIVITY','yes', ...
                    'TOL',1e-4, 'AL_iters',2000,'verbose','yes');
timeSunsal = toc;


%% Sunsal total variation

lambda = 1e-3;
lambda_TV = 1e-3;

tic
[X_hat_tv,res,rmse] = sunsal_tv(A,Y,'MU',0.05,'POSITIVITY','yes','ADDONE','no', ...
                              'LAMBDA_1',lambda,'LAMBDA_TV', lambda_TV, 'TV_TYPE','niso',...
                              'IM_SIZE',[nl,nc],'AL_ITERS',200,  'VERBOSE','yes');
timeTV = toc;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constrained least squares using l2-l1-swSp solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda_swsp = 0.07e-3; 0.05e-3; 0.002;

XT = zeros(size(X_hat_tv));

tic
[X_hat_l21LC,res2,rmse_ni] = sunsal_tv_lw_sp(A,Y,'MU',0.5,'POSITIVITY','yes','ADDONE','no', ...
                               'LAMBDA_1',lambda_swsp,'IM_SIZE',[nl,nc],'AL_ITERS',5,'TRUE_X', XT, 'VERBOSE','yes');  
time_swsp = toc;










%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print resulting time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% timeKmeans = -1;

fprintf('\n\n EXECUTION TIME \n\n')

fprintf('TIME-l1 = %2.3f\n', timeSunsal)
fprintf('TIME-TV = %2.3f\n', timeTV)
fprintf('TIME_SPx-Mscale = %2.3f\n', timeMscale)
fprintf('TIME-kmeans = %2.3f\n', timeKmeans)
fprintf('TIME-BPT = %2.3f\n', timeBPT)
fprintf('TIME-SWSP = %2.3f\n\n', time_swsp)  
        

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot abundances per pixel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(131)
imagesc(X_hat_l1)
axis off
title('SUnSAL')

subplot(132)
imagesc(X_hat_tv)
axis off
title('SUnSAL-TV')

subplot(133)
imagesc(X_hat_l1_spreg)
axis off
title('Mscale-SprPx')




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot abundance maps for the desired materials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Alunite-s in the lib: 231, 232, 256, 257, 420, 456
% Chalcedony-ies in the lib: 297
% Muscovite-s  in the lib: 46, 47, 93, 94, 95, 108, 153, 156, 172, 177, 203, 212, 253
% Buddingtonite-s in the lib: 336, 337
% Dickite-s in the lib: 185, 236
% Calcite-s in the lib: 34, 35, 81


% material_idx = 420; % alunite
% material_idx = 336; % Buddingtonite
% material_idx = 203; % muscovite
% material_idx = 297; % Chalcedony


%%

X_hat_l1_im = reshape(X_hat_l1', nl,nc,n); 
X_hat_tv_im = reshape(X_hat_tv', nl,nc,n); 
X_hat_l1_spreg_im = reshape(X_hat_l1_spreg', nl,nc,n);
X_hat_l1_BPT_im = reshape(X_hat_l1_BPT', nl,nc,n); 
X_hat_l21LC_im = reshape(X_hat_l21LC', nl,nc,n); 


material_idx = [420 336 297];
range_mat = [0 0.38; ... % Alunite
             0 0.30; ... % Buddingtonite
             0 0.65];    % Chalcedony
name_mat = {'Alunite','Buddingtonite','Chalcedony'};

count = 0;
for idx=material_idx
    count = count + 1;
    
    figure
%     [ha, pos] = tight_subplot(1,4,[.025 .025],[.05 .05],[.05 .15]);
    [ha, pos] = tight_subplot(1,4,[.025 .015],[.05 .05],[.05 .15]);
    
    ii = 1;

    axes(ha(ii)); ii = ii+1;
    imagesc(X_hat_tv_im(:,:,idx), range_mat(count,:))
    axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])

    axes(ha(ii)); ii = ii+1;
    imagesc(X_hat_l21LC_im(:,:,idx), range_mat(count,:))
    axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])

    axes(ha(ii)); ii = ii+1;
    imagesc(X_hat_l1_BPT_im(:,:,idx), range_mat(count,:))
    axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])

    axes(ha(ii)); ii = ii+1;
    imagesc(X_hat_l1_spreg_im(:,:,idx), range_mat(count,:))
    axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
    originalSize2 = get(gca, 'Position');

    h=colorbar; 
    set(ha(4), 'Position', originalSize2);
    % set(h,'fontsize',5);
    set(h,'fontsize',8);

    colormap jet
    print(strcat('examples/realImg/estim_abundances_',name_mat{count},'_tght'),'-dpdf')


end





