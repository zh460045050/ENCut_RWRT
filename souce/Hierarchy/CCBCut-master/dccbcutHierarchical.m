function [I,content,cost,t] = dccbcutHierarchical(varargin)
% dccbcutHierarchical: computes hierarchical CCBCut with downsampling of 
%   weighted adjacency matrix
% usage: I = dccbcutHierarchical(W,imSize,numDownsamp,k,options);
%    or: [I,content,cost,t] = dccbcutHierarchical(...);
% arguments:
%   W - adjacency matrix, assumed to be constructed from an mxn image so
%       that the rows of W correspond to pixels in the image in 
%       column-major order 
%   imSize - 1x2 array containing number of rows and columns of image that
%       was used to construct W
%   numDownsamp - number of downsampling operations to perform (default
%       numDownsamp = 2)
%   k - desired number of subgraphs (default k = 2)
%   options - argument created with the ccbcutSet function. See ccbcutSet 
%       for details. 
%
%   I - array containing indices into resulting subgraphs for each vertex
%   content (kx2) - array containing information about partitions:
%       content(:,1) has number of vertices in each partition
%       content(:,2) has total degree of each partition
%   cost - array (or structure) containing cost at each iteration
%   t - array (or structure) containing total computational effort for
%       1-spectral clustering step

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is a modified version of the DNCuts code for fast eigenvector
% computation in Normalized Cuts Segmentation, from the Berkeley group,
% based on the paper: "Multiscale Combinatorial Grouping," P. Arbelaez et
% al., CVPR 2014.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% author: Nathan D. Cahill
% email: nathan.cahill@rit.edu
% date: 21 February 2018

% parse input arguments
[W,imSize,numDownsamp,numClusters,opts] = parseInputs(varargin{:});

% initialize variables for downsampling affinity matrix
Wd = W;
imSizeD = imSize;
Bs = cell(numDownsamp,1);

% downsample affinity matrix (and Y0, if necessary)
for i = 1:numDownsamp
  
  % Create a binary array of the pixels that will remain after decimating
  
  % every other row and column
  [r,c] = ind2sub(imSizeD, 1:size(Wd,1));
  do_keep = (mod(r, 2) == 0) & (mod(c, 2) == 0);
  
  % Downsample the affinity matrix
  Ws = Wd(:,do_keep)';
  
  % Normalize the downsampled affinity matrix
  d = (sum(Ws,1) + eps);
  B = bsxfun(@rdivide, Ws, d)';
  
  % "Square" the affinity matrix, while downsampling
  Wd = Ws*B;

  imSizeD = floor(imSizeD / 2);
  
  % Hold onto the normalized affinity matrix for bookkeeping
  Bs{i} = B;  

end

% force symmetry of Wd
Wd = (Wd + Wd')/2;

% zero out main diagonal elements of Wd to avoid self-loops
for j = 1:size(Wd,1)
    Wd(j,j) = 0;
end

% perform multipartitioning on downsampled adjacency matrix
tic; [Id,cuts] = computeMultiPartitioning_ccbcut(Wd,opts,numClusters); t = toc;
Id = Id(:,end);
cost = cuts(end);

% Upsample the resulting clustering
%I = imresize(reshape(Id,imSizeD([2 1])),imSize([2 1]),'nearest')';
I = imresize(reshape(Id,imSizeD),imSize,'nearest');

% Compute content of each partition
content = zeros(numClusters,2);
content(:,1) = histc(I(:),1:numClusters);
degVals = full(sum(W,2));
for j = 1:numClusters
    content(j,2) = sum(degVals(I==j));
end

end

function [W,imSize,numDownsamp,numClusters,opts] = parseInputs(varargin)
% parse input arguments

% check number of inputs
narginchk(2,inf);

% get/check W
W = varargin{1};
if ~isreal(W) || ~isequal(ndims(W),2) || ~isequal(size(W,1),size(W,2))
    error([mfilename,':InvalidW'],...
        'W must be real square weight matrix.');
end

% get/check imSize
imSize = varargin{2};
if ~isnumeric(imSize) || ~isequal(numel(imSize),2) || ~isequal(prod(imSize),size(W,1))
    error([mfilename,':InvalidK'],...
        'imSize must be 1x2 array whose product is the number of rows in W.');
end

% get/check numDownsamp
if nargin<3
    numDownsamp = [];
else
    numDownsamp = varargin{3};
end
if isempty(numDownsamp)
    numDownsamp = 2;
end
if ~isscalar(numDownsamp) || ~isequal(numDownsamp,round(numDownsamp)) || numDownsamp<1
    error([mfilename,':InvalidNumDownsamp'],...
        'numDownsamp must be positive integer.');
end

% get/check numClusters
if nargin<4
    numClusters = [];
else
    numClusters = varargin{4};
end
if isempty(numClusters)
    numClusters = 2;
end
if ~isscalar(numClusters) || ~isequal(numClusters,round(numClusters)) || numClusters<2
    error([mfilename,':InvalidNumClusters'],...
        'numClusters must be positive integer >= 2.');
end

% get/check opts
if nargin<5
    opts = [];
else
    opts = varargin{5};
end
if isempty(opts)
    opts = ccbcutSet('Display','final');
end

end

