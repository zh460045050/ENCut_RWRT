function [I,Y,content,cost,t,costOrthCon,gamma,e] = dccbcut(varargin)
% dccbcut: computes CCBCut with downsampling of weighted adjacency matrix
% usage: I = dccbcut(W,imSize);
%    or: I = dccbcut(W,imSize,numDownsamp);
%    or: I = dccbcut(W,imSize,numDownsamp,k);
%    or: I = dccbcut(W,imSize,numDownsamp,k,Y0);
%    or: I = dccbcut(W,imSize,numDownsamp,k,Y0,options);
%    or: [I,Y] = dccbcut(...);
%    or: [I,Y,content] = dccbcut(...);
%    or: [I,Y,content,cost] = dccbcut(...);
%    or: [I,Y,content,cost,t] = dccbcut(...);
%    or: [I,Y,content,cost,t,costOrthCon] = dccbcut(...);
%    or: [I,Y,content,cost,t,costOrthCon,gamma] = dccbcut(...);
%    or: [I,Y,content,cost,t,costOrthCon,gamma,e] = dccbcut(...);
% 
% arguments:
%   W - adjacency matrix, assumed to be constructed from an mxn image so
%       that the rows of W correspond to pixels in the image in 
%       column-major order 
%   imSize - 1x2 array containing number of rows and columns of image that
%       was used to construct W
%   numDownsamp - number of downsampling operations to perform (default
%       numDownsamp = 2)
%   k - desired number of subgraphs (default k = 2)
%   Y0 - size(W,1) x (k-1) array containing starting guess for embedding
%       coordinates. (Alternatively, Y0 can be a cell array, the last
%       element of which contains a starting guess at the coarsest level.
%       No error checking will be done to make sure sizes are consistent in
%       this case, so use this option with caution!)
%   options - argument created with the ccbcutSet function. See ccbcutSet for
%       details. 
%
%   I - array containing indices into resulting subgraphs for each vertex
%   Y - cell array containing (k-1)-dimensional embedding coordinates for 
%       each vertex at each level of coarsening
%   content (kx2) - array containing information about partitions:
%       content(:,1) has number of vertices in each partition
%       content(:,2) has total degree of each partition
%   cost - array (or structure) containing cost at each iteration
%   t - array (or structure) containing time spent in each iteration
%   costOrthCon - array containing orthogonality constraint costs at each
%       iteration (only for Bregman algorithm - returns empty for IRRQ)
%   gamma - optimal weights from IRRQ solution (only for IRRQ algorithm -
%       returns empty for Bregman)
%   e - epsilon value at final iteration of IRRQ (only for IRRQ algorithm -
%       returns empty for Bregman)

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
[W,imSize,numDownsamp,ccbcutArgs,Y0flag,Y0coarseGridFlag,defaultopts] = parseInputs(varargin{:});

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
  
  % Downsample Y0 if provided on finest level
  if Y0flag && ~Y0coarseGridFlag
      ccbcutArgs{2} = ccbcutArgs{2}(do_keep,:);
  end
  
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

% if cell array was passed for Y0, strip out initialization and pass
% directly to ccbcut as numeric array
if Y0flag && Y0coarseGridFlag
    ccbcutArgs{2} = ccbcutArgs{2}{numDownsamp+1};
end

% determine which type of rounding to use, and set ccbcut options to do no
% rounding
if numel(ccbcutArgs)>=3
    opts = ccbcutArgs{3};
else
    opts = defaultopts;    
end
roundingType = ccbcutGet(opts,'RoundingType',defaultopts,'fast');
ccbcutArgs{3} = ccbcutSet(opts,'RoundingType','none');
switch roundingType
    case 'kmeans'
        numStartPoints = ccbcutGet(opts,'KmeansNumStartPoints',defaultopts,'fast');
        kmeansDistance = ccbcutGet(opts,'KmeansDistance',defaultopts,'fast');
    case 'graclus'
        if ~isunix
            error('graclus rounding only supported on UNIX.');
        end
end

% perform ccbcuts on downsampled adjacency matrix
Y = cell(numDownsamp+1,1);
[~,Y{numDownsamp+1},~,cost,t,costOrthCon,gamma,e] = ccbcut(Wd, ccbcutArgs{:});

% Upsample the resulting embedding
for i = numDownsamp:-1:1
  Y{i} = Bs{i} * Y{i+1};
end

% whiten the eigenvectors, as they can get scaled weirdly during upsampling
%Y = whiten(Y, 1, 0);

% round embedding to form partitioning
I = []; content = [];
switch roundingType
    case 'kmeans'
        
        if numStartPoints>0
            k = size(Y{1},2)+1;
            fprintf(1,'Performing K-means rounding...\n');
            I = kmeansRounding(Y{1},numStartPoints,k,kmeansDistance);
            
            % compute content of each partition
            content = zeros(k,2);
            content(:,1) = histc(I,1:k);
            degVals = full(sum(W,2));
            for j = 1:k
                content(j,2) = sum(degVals(I==j));
            end
        end

    case 'graclus'
        
        % find indices of nonzero weights
        [i,j,wVals] = find(W);
        
        % compute gamma weights
        tau = ccbcutGet(opts,'Tau',defaultopts,'fast');
        gamma = ones(numel(i),1);
        if tau<2
            for m = 1:numel(i)
                Yi = Y{1}(i(m),:);
                Yj = Y{1}(j(m),:);
                gamma(m) = (wVals(m).*sum((Yi-Yj).^2) + e.^2)^((tau-2)/2);
            end
        elseif tau>2
            for m = 1:numel(i)
                Yi = Y{1}(i(m),:);
                Yj = Y{1}(j(m),:);
                gamma(m) = (wVals(m).*sum((Yi-Yj).^2))^((tau-2)/2);
            end
        end
        gamma = gamma./mean(gamma);
        
        % construct diagonal normalization matrix Pi
        D = diag(sum(W,2));
        switch ccbcutGet(opts,'BalanceType',defaultopts,'fast')
            case 'normalized'
                Pi = D;
            case 'ratio'
                Pi = speye(size(W));
        end

        % construct k
        k = size(Y{1},2)+1;
        
        % perform rounding
        fprintf(1,'Performing Graclus rounding...\n');
        I = graclusRounding(W,gamma,Pi,k);
        
        % compute content of each partition
        content = zeros(k,2);
        content(:,1) = histc(I,1:k);
        degVals = full(sum(W,2));
        for j = 1:k
            content(j,2) = sum(degVals(I==j));
        end
        
end

end

function [W,imSize,numDownsamp,ccbcutArgs,Y0flag,Y0coarseGridFlag,defaultopts] = parseInputs(varargin)
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

% the rest of the arguments will passed directly to ccbcut and tested there
ccbcutArgs = varargin(4:end);

% if Y0 was provided, will need to downsample it
Y0flag = false; Y0coarseGridFlag = false;
if numel(ccbcutArgs)>=2 && ~isempty(ccbcutArgs{2})
    Y0flag = true;
    if isa(ccbcutArgs{2},'cell') && isequal(numel(ccbcutArgs{2}),numDownsamp+1)
        Y0coarseGridFlag = true;
    end
end

defaultopts = struct('Display','final','Algorithm','irrq','Tau',1,...
    'BalanceType','normalized','RoundingType','kmeans',...
    'KmeansNumStartPoints',10,'KmeansDistance','sqeuclidean',...
    'MaxSOCIters',100,'MaxL1ProbIters',100,'TolBregSOCX',1e-3,...
    'TolBregSOCF',1e-3,'TolBregL1ProbX',1e-3,'TolBregL1ProbF',1e-3,...
    'BregLambda',1000,'R',100,'MinIRRQIters',5,'MaxIRRQIters',100,'TolIRRQX',1e-3,'TolIRRQF',1e-3,...
    'TolIRRQL1',1e-3,'TolIRRQL1Curv',1e-1,'Epsilon',1,'KRatio',0.5,'IRRQConstraint','none',...
    'IRRQLambda',1,'IRRQq',1,'IRRQAlpha',0.5,'EpsilonDenom',1,...
    'KRatioDenom',0.5); 

end

