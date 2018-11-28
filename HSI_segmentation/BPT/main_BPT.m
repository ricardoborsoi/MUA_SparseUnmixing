%% Main for BPT construction and processing
%% Initial segmentation

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

%% build data structures

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


%% build BPT

T = buildBPT(DAT);

%% Process/prune BPT

% BPT processing can be done in different way, and must be specifically
% tuned to the desired goal.

% A simple way to process the BPT is to extract a given number of regions N,
% which will constitute a partition featuring the N most dissimilar regions
% created during the construction of the BPT.

prunedtree = pruneBPTnbregions(T,50); % N=10

% A more elaborated processing method is to design some energetic function
% and to extract the partition minimizing the overall energy (typically a
% Mumford-Shah energetic functional).
% The region-wise energy has to be composed of a goodness-of-fit term and a
% regularization term, with a trade-off parameter 'lambda'. 
% If the energy of the partition is the sum of the energies of the regions 
% composing the partition, then it is possible to extract at once all
% optimal partitions with respect to the value of 'lambda'
% It is then possible to extract an optimal partiion either by specifiying
% a value for the 'lambda' regularization parameter (option "unconstrained"
% in the pruneBPTminimizeEnergy), by specifying a bound for the total
% energy value of the partition (option "constrained") or by specifying a
% certain number of regions (option "region"), the returned optimal
% partition being as close as possible from the given number of regions

% % % prunedtree = getBPTenergy(T,DAT.data,DAT.initsegmap,@GofF_MS,@Regu_boundary);
% % % prunedtree = getpersistentintervals(prunedtree,'display');
% % % prunedtree = pruneBPTminimizeEnergy(prunedtree,250,'regions');

%% Display partition

segregions = retrievesegmentation(prunedtree,initsegmap,'regions',imrgb);
drawborders(imrgb,segregions,'red');
imrgbmean = displaysegmentationfalsecolors(segregions,imrgb);
figure, imshow(imrgbmean)
