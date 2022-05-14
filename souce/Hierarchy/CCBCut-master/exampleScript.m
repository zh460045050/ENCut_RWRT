%% update search path to include code directory and subdirectory

%%%%%%%%%%%%%%%%%%
% If your current working directory is the 'code' directory, this line is
% fine. If not, change this line to indicate the absolute path of the 
% 'code' directory.
codeDir = pwd;  
%%%%%%%%%%%%%%%%%%

addpath(genpath(codeDir));

%% load image to segment
img = imread('2018.jpg');

%% convert to L*a*b*
cform = makecform('srgb2lab');
imgLab = im2double(applycform(img,cform));

%% display image
figure; imshow(img); title('Original Image');

%% construct affinity matrix

% pixel radius for including nonzero weights
r = 10;

% construct indices of neighbors
[x,y] = ndgrid(-r:r,-r:r);
d = (x.^2 + y.^2) <= r.^2;
d(r+1,r+1) = false;
neighbors = [x(d) y(d)];

% construct affinity matrix
fprintf(1,'Constructing affinity matrix...\n\t(If r=10, this may take ~2 minutes...)\n');
tic; W = constructWeightMatrix(imgLab,neighbors,90); toc
W = (W+W')/2;

%% desired number of clusters
clusterNumber = 25;

%% desired value of tau
tau = 1;

%% set up parameters for optimization

% parameters for CCBCuts
kRatio = 0.0001;
minItersNCut = 5;
minItersRCut = 2;
maxIters = 50;
tolVal = eps;
tolL1 = 1e-2;
tolCurvPhaseI = 0.08;
tolCurvPhaseII = 1;
numDownsamp = 3; % 2 used in paper; use 3 or 4 for faster but less accurate results
numStartPoints = 1;
imSize = [size(img,1) size(img,2)];

% options structure for CCNCuts
optsCCNCut = ccbcutSet('Display','iter',...
    'Algorithm','irrq',...
    'Tau',tau,...
    'BalanceType','normalized',...
    'RoundingType','kmeans',...
    'KmeansNumStartPoints',numStartPoints,...
    'MinIRRQIters',minItersNCut,...
    'MaxIRRQIters',maxIters,...
    'TolIRRQX',tolVal,...
    'TolIRRQF',tolVal,...
    'TolIRRQL1',tolL1,...
    'TolIRRQL1Curv',tolCurvPhaseI,...
    'Epsilon',1,...
    'KRatio',kRatio,...
    'IRRQConstraint','hard');

% options for CCRCuts
optsCCRCut = ccbcutSet('Display','iter',...
    'Algorithm','irrq',...
    'Tau',tau,...
    'BalanceType','ratio',...
    'RoundingType','kmeans',...
    'KmeansNumStartPoints',numStartPoints,...
    'MinIRRQIters',minItersRCut,...
    'MaxIRRQIters',maxIters,...
    'TolIRRQX',tolVal,...
    'TolIRRQF',tolVal,...
    'TolIRRQL1',tolL1,...
    'TolIRRQL1Curv',tolCurvPhaseI,...
    'Epsilon',1,...
    'KRatio',kRatio,...
    'IRRQConstraint','hard');

%% perform Multiway CCNCuts
fprintf(1,'Multiway CCNCut with IRRQ, tau = %g...\n',tau);
[idx,Y] = dccbcut(W,imSize,numDownsamp,clusterNumber,[],optsCCNCut);
imgSegLabels = uint8(reshape(idx,imSize));
imgSegCCNCut = colorSegmentedImage(img,imgSegLabels);

% display embedding
figure;
for j = 1:clusterNumber-1
    subplot(4,6,j);
    imagesc(reshape(Y{1}(:,j),imSize));
    axis image; set(gca,'xtick',[],'ytick',[]);
end 
set(gcf,'name',sprintf('Embedding Components from Multiway CCNCuts algorithm, tau = %g',tau));
    
% display segmentation result
figure; 
subplot(1,2,1); imagesc(imgSegLabels); axis image; 
set(gca,'xtick',[],'ytick',[]);
title('Labels');
subplot(1,2,2); imshow(imgSegCCNCut);
title('Colorized');
set(gcf,'name',sprintf('Segmentation using Multiway CCNCuts algorithm, tau = %g',tau));

%% perform Multiway CCRCuts
fprintf(1,'Multiway CCRCut with IRRQ, tau = %g...\n',tau);
[idx,Y] = dccbcut(W,imSize,numDownsamp,clusterNumber,[],optsCCRCut);
imgSegLabels = uint8(reshape(idx,imSize));
imgSegCCRCut = colorSegmentedImage(img,imgSegLabels);

% display embedding
figure;
for j = 1:clusterNumber-1
    subplot(4,6,j);
    imagesc(reshape(Y{1}(:,j),imSize));
    axis image; set(gca,'xtick',[],'ytick',[]);
end 
set(gcf,'name',sprintf('Embedding Components from Multiway CCRCuts algorithm, tau = %g',tau));

% display segmentation result
figure; 
subplot(1,2,1); imagesc(imgSegLabels); axis image; 
set(gca,'xtick',[],'ytick',[]);
title('Labels');
subplot(1,2,2); imshow(imgSegCCRCut);
title('Colorized');
set(gcf,'name',sprintf('Multiway Segmentation using CCRCuts algorithm, tau = %g',tau));

%% perform Hierarchical 2-way CCNCuts
fprintf(1,'Hierarchical 2-way CCNCut with IRRQ, tau = %g...\n',tau);
idx = dccbcutHierarchical(W,imSize,numDownsamp,clusterNumber,optsCCNCut);
imgSegLabels = uint8(reshape(idx,imSize));
imgSegCCNCut = colorSegmentedImage(img,imgSegLabels);

% display result
figure; 
subplot(1,2,1); imagesc(imgSegLabels); axis image; 
set(gca,'xtick',[],'ytick',[]);
title('Labels');
subplot(1,2,2); imshow(imgSegCCNCut);
title('Colorized');
set(gcf,'name',sprintf('Hierarchical 2-Way Segmentation using CCNCuts algorithm, tau = %g',tau));

%% perform Hierarchical 2-way CCRCuts
fprintf(1,'Hierarchical 2-way CCRCut with IRRQ, tau = %g...\n',tau);
idx = dccbcutHierarchical(W,imSize,numDownsamp,clusterNumber,optsCCRCut);
imgSegLabels = uint8(reshape(idx,imSize));
imgSegCCNCut = colorSegmentedImage(img,imgSegLabels);

% display result
figure; 
subplot(1,2,1); imagesc(imgSegLabels); axis image; 
set(gca,'xtick',[],'ytick',[]);
title('Labels');
subplot(1,2,2); imshow(imgSegCCNCut);
title('Colorized');
set(gcf,'name',sprintf('Hierarchical 2-Way Segmentation using CCRCuts algorithm, tau = %g',tau));
