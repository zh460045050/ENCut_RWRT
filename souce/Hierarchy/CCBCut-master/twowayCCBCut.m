

img = imread('imgs/100007.jpg');
img = imresize(img, 0.5);
tau = 1;
%% convert to L*a*b*
cform = makecform('srgb2lab');
imgLab = im2double(applycform(img,cform));

% r = 2;
% [x,y] = ndgrid(-r:r,-r:r);
% d = (x.^2 + y.^2) <= r.^2;
% d(r+1,r+1) = false;
% neighbors = [x(d) y(d)];
% W = constructWeightMatrix(imgLab,neighbors,90);
% W = (W+W')/2;



sigma = 50;
[X, Y, Z] = size(imgLab);
N = X * Y;
imgVals = reshape(imgLab,N,Z); clear Lab;
[~, edges] = lattice(X,Y,1);
weights = makeweights_old(edges,imgVals,sigma);
W = adjacency(edges,weights,N); clear edges weights;


numDownsamp = 1; % 2 used in paper; use 3 or 4 for faster but less accurate results
k = 3;
imSize = [size(img,1) size(img,2)];

optsCCNCut = ccbcutSet('Display','iter',...
    'Algorithm','irrq',...
    'Tau',1,...
    'BalanceType','normalized',...
    'RoundingType','kmeans',...
    'KmeansNumStartPoints',1,...
    'MinIRRQIters',5,...
    'MaxIRRQIters',50,...
    'TolIRRQX',eps,...
    'TolIRRQF',eps,...
    'TolIRRQL1',1e-2,...
    'TolIRRQL1Curv',0.08,...
    'Epsilon',1,...
    'KRatio',0.0001,...
    'IRRQConstraint','hard');

fprintf(1,'Hierarchical 2-way CCNCut with IRRQ, tau = %g...\n',tau);
idx = dccbcutHierarchical(W,imSize,numDownsamp,k,optsCCNCut);
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


result = show_result(img, imgSegLabels);
figure();
imshow(result);