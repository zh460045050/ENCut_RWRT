function [labels,Y] = ccbcutSegment(varargin)
% ccbcutSegment: perform image segmentation using CCBCut embeddings
% usage: labels = ccbcutSegment(img,seeds,neighbors);
%    or: [labels,Y] = ccbcutSegment(img,seeds,neighbors,labels0);
%
% arguments:
%   img (mxnx3) - color image
%   seeds (mxn) - array containing indices of user-provided scribbles.
%       Scribble pixels should be labeled from 0, 1, ..., k-1. Non-scribble
%       pixels should be labeled -1.
%   neighbors (rx2) - array containing integer offsets from (0,0)
%       indicating positions of neighbors to be used in graph construction
%       Default neighbors = [-1 0;0 -1;1 0;0 1];
%   labels0 (mxn) - initial guesses for class labels. If not provided, 
%       GMM will be used to initialize.
%
%   labels (mxn) - class labels
%   Y (mxnxk) - embedding coordinates
%

% author: Nathan D. Cahill
% email: nathan.cahill@rit.edu
% date: 21 February 2018

% parse input arguments
[img,~,neighbors,labels0,k] = parseInputs(varargin{:});

% convert image to double
img = im2double(img);
imgVec = reshape(img,[size(img,1)*size(img,2),3]);

% construct weighted adjacency matrix from image
W = constructWeightMatrix(img,neighbors);

% initialize labels with GMM if none (or all zeros) provided
if ~any(labels0(:))
    
    gmmOpts = statset('Display','Off','MaxIter',200);
    regVal = 0.0001;
    
    GMModel = fitgmdist(imgVec,k,...
        'Options',gmmOpts,'RegularizationValue',regVal);
    
    [idx,~] = cluster(GMModel,imgVec);
    
    labels0 = reshape(idx,[size(img,1),size(img,2)]);

else % add 1 to labels so they range from 1 ... k
    
    labels0 = labels0 + 1;
    
end

% use initial labels to generate initial estimate of embedding
Y0 = sparse((1:size(imgVec,1))',labels0,1,size(imgVec,1),k);
Y0 = full(Y0./repmat(sqrt(sum(Y0.^2)),[size(imgVec,1),1]));

% first stage of piecewise flat embedding
opts1 = ccbcutSet('Display','notify',...
    'MaxSOCIters',10,...
    'MaxL1ProbIters',5,...
    'TolSOC',1e-6,...
    'TolL1Prob',1e-6,...
    'Lambda',10000,...
    'R',100);
Y1 = ccbcut(W,k,Y0,opts1);

% second stage of piecewise flat embedding
opts2 = ccbcutSet('Display','iter',...
    'MaxSOCIters',1,...
    'MaxL1ProbIters',100,...
    'TolSOC',1e-4,...
    'TolL1Prob',1e-6,...
    'Lambda',10000,...
    'R',10);
Y2 = ccbcut(W,k,Y1,opts2);

% k-means clustering on result
labels2 = bestPermutation(labels0(:),kmeans(Y2,k));

% reshape results 
Y = reshape(Y2,[size(img,1),size(img,2),k]);
labels = reshape(labels2,[size(img,1),size(img,2)]);

% subtract 1 from labels so that they range from 0 ... k-1
labels = labels - 1;

end

function [img,seeds,neighbors,labels0,numClusters] = parseInputs(varargin)
% parse input arguments

% check number of inputs
narginchk(2,4);

% get/check img
img = im2double(varargin{1});

% get/check seeds
seeds = varargin{2};
if ~isequal(size(seeds),[size(img,1),size(img,2)]) || ~isequal(seeds,round(seeds)) || ~isequal(unique(seeds(:)),(-1:max(seeds(:)))')
    error([mfilename,':InvalidSeeds'],...
        'seeds must be 2-D array of integers from -1 to k-1, same number of rows and columns as img.');
end

% initialize number of clusters based on seeds
numClusters = max(seeds(:))+1;

% get/check neighbors
if nargin<3
    neighbors = [];
else
    neighbors = varargin{3};
end
if isempty(neighbors)
    neighbors = [-1 0;0 -1;1 0;0 1];
end
if ~isequal(ndims(neighbors),2) || ~isequal(size(neighbors,2),2) || ~isequal(neighbors,round(neighbors))
    error([mfilename,':InvalidNeighbors'],...
        'Neighbors must be rx2 array of integers.');
end

% get/check labels0
if nargin<4
    labels0 = [];
else
    labels0 = varargin{4};
end
if isempty(labels0)
    labels0 = zeros(size(img,1),size(img,2));
end
if ~isequal(size(labels0),[size(img,1),size(img,2)]) || any(labels0(:)<0) || ~isequal(labels0,round(labels0))
    error([mfilename,':InvalidLabels0'],...
        'labels0 must be 2-D array of nonnegative integers, same number of rows and columns as img.');
end

% update numClusters if extra labels have been provided in labels0
if max(labels0(:)) > numClusters - 1
    numClusters = max(labels0(:)) + 1;
end

end

function W = constructWeightMatrix(img,neighbors)
% construct weight matrix for image

% number of rows and columns in image
r = size(img,1); c = size(img,2);

% pad image by one-pixel boundary
padSize = max(abs(neighbors));
imgPad = padarray(img,padSize,'replicate','both');
imgRows = (1:r)+padSize(1);
imgCols = (1:c)+padSize(2);
maskInd = find(padarray(true(r,c),padSize,false));
[maskI,maskJ] = ind2sub([r+2*padSize(1),c+2*padSize(2)],maskInd);

% loop through neighbors, constructing sparse matrix of weights for each
% specific neighbor individually
numNeighbors = size(neighbors,1);
dMax = zeros(numNeighbors,1);
A = cell(numNeighbors,1);
beta = 90;
numNonZeros = zeros(numNeighbors,1);

% minimum weight for each neighbor
dMin = 1e-4;

for i = 1:numNeighbors
    d = sqrt(sum((imgPad(imgRows+neighbors(i,1),imgCols+neighbors(i,2),:)-img).^2,3));
    d = max(d,dMin);
    dMax(i) = max(d(:));
    
    % compute linear indices of neighbors
    neighborsInd = sub2ind([r+2*padSize(1),c+2*padSize(2)],...
        maskI + neighbors(i,1),maskJ + neighbors(i,2));
    
    % construct sparse matrix containing these neighbors
    A{i} = sparse((1:length(maskInd))',neighborsInd,d,...
        length(maskInd),(r+2*padSize(1))*(c+2*padSize(2)));

    % remove columns corresponding to padded elements
    A{i} = A{i}(:,maskInd);
    
    % determine number of nonzero entries
    numNonZeros(i) = nnz(A{i});
    
end

% now normalize weights, and compute exp(-beta*weight) for all nonzero
% weights
dMaxAll = max(dMax);
for i = 1:numNeighbors
    
    % normalize weights
    A{i} = A{i}./dMaxAll;
    
    % compute exp(-beta*weight)
    [rowInd,colInd,val] = find(A{i});
    [numRows,numCols] = size(A{i});
    A{i} = sparse(rowInd,colInd,exp(-beta.*val),numRows,numCols);
    
end

% add each of the W{i}'s together to form final weight matrix
W = spalloc(r*c,r*c,sum(numNonZeros));
for i = 1:numNeighbors
    W = W + A{i};
end

end