function [I,Y,content,cost,t,costOrthCon,gamma,e,numStartPoints,kmeansDistance] = ccbcut(varargin)
% ccbcut: computes CCBCut given weighted adjacency matrix
% usage: I = ccbcut(W,k);
%    or: I = ccbcut(W,k,Y0);
%    or: I = ccbcut(W,k,Y0,options);
%    or: [I,Y] = ccbcut(...);
%    or: [I,Y,content] = ccbcut(...);
%    or: [I,Y,content,cost] = ccbcut(...);
%    or: [I,Y,content,cost,t] = ccbcut(...);
%    or: [I,Y,content,cost,t,costOrthCon] = ccbcut(...);
%    or: [I,Y,content,cost,t,costOrthCon,gamma] = ccbcut(...);
%    or: [I,Y,content,cost,t,costOrthCon,gamma,e] = ccbcut(...);
% 
% arguments:
%   W - adjacency matrix
%   k - desired number of subgraphs (default k = 2)
%   Y0 - size(W,1) x (k-1) array containing starting guess for embedding
%       coordinates
%   options - argument created with the ccbcutSet function. See ccbcutSet for
%       details. 
%
%   I - array containing indices into resulting subgraphs for each vertex
%   Y - array containing (k-1)-dimensional embedding coordinates for each
%       vertex
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

% author: Nathan D. Cahill
% email: nathan.cahill@rit.edu
% date: 21 February 2018

% parse input arguments
[W,k,Y,options,defaultopts] = parseInputs(varargin{:});

% initialize outputs that might not be set
cost = []; t = []; costOrthCon = []; gamma = []; e = [];

% get option settings
displayFlag = ccbcutGet(options,'Display',defaultopts,'fast');
balanceType = ccbcutGet(options,'BalanceType',defaultopts,'fast');
roundingType = ccbcutGet(options,'RoundingType',defaultopts,'kmeans');
tau = ccbcutGet(options,'Tau',defaultopts,'fast');
alg = ccbcutGet(options,'Algorithm',defaultopts,'fast');
numStartPoints = ccbcutGet(options,'KmeansNumStartPoints',defaultopts,'fast');
kmeansDistance = ccbcutGet(options,'KmeansDistance',defaultopts,'fast');

% error out for unsupported combinations of parameters
switch alg
    case 'irrq'
%         if tau>2
%             error('IRRQ algorithm not currently supported for tau>2.');
%         end
    case 'bregman'
        if ~isequal(tau,2) || ~isequal(tau,1)
            error('Bregman algorithm not currently supported for tau not equal to 1.');
        end
end
switch roundingType
    case 'graclus'
        if ~isunix
            error('Graclus rounding only currently supported in UNIX versions of MATLAB.');
        end
end

% construct diagonal normalization matrix 
D = diag(sum(W,2));
switch balanceType
    case 'normalized'
        Pi = D;
    case 'ratio'
        Pi = speye(size(W));
end

% initial estimate of Y 
if isempty(Y) % compute intial estimate using Laplacian Eigenmaps
    
    % display if necessary
    switch displayFlag
        case {'iter','notify'}
            fprintf(1,'Computing initial estimate of embedding using Laplacian Eigenmaps...\n');
    end
    
    % construct Laplacian matrix
    L = D - W;
    try
        [Y,~] = lobpcg(randn(size(W,1),k),L,Pi,[],1e-5,100*size(W,1));
    catch
        [Y,~] = eigs(L,Pi,k,'SA');
    end

    % discard eigenvector corresponding to zero eigenvalue
    Y = Y(:,2:k);

    % ensure unit Pi-norm
    for m = 1:(k-1)
        Y(:,m) = Y(:,m)./sqrt(Y(:,m)'*(Pi*Y(:,m)));
    end
    
end

% initialize gamma weights
gamma = ones(nnz(W),1);

% update estimate if tau is not 2
if ~isequal(tau,2)
    
    % perform either IRRQ or Bregman algorithm
    switch alg

        case 'irrq'

            % display if necessary
            switch displayFlag
                case {'iter','notify'}
                    fprintf(1,'Performing IRRQ algorithm to compute CCBCut embedding...\n');
            end
            
            [Y,cost,t,gamma,e] = ccbcutIRRQ(W,k-1,Y,options,defaultopts);

        case 'bregman'

            % display if necessary
            switch displayFlag
                case {'iter','notify'}
                    fprintf(1,'Performing Bregman algorithm to compute CCBCut embedding...\n');
            end
            
            [Y,cost,t,costOrthCon] = ccbcutBregman(W,k-1,Y,options,defaultopts);

    end

else
   
    % display if necessary
    switch displayFlag
        case {'iter','notify','final'}
            fprintf(1,'No iterations necessary; initial estimate is exact for tau = 2...\n');
    end
        
end

% display if necessary
switch displayFlag
    case {'iter','notify','final'}
        switch roundingType
            case 'kmeans'
                fprintf(1,'Rounding embedding using Kmeans clustering...\n');
            case 'graclus'
                fprintf(1,'Rounding embedding using Graclus...\n');

        end
end

% round embedding to form partitioning
I = []; content = [];
switch roundingType
    case 'kmeans'
        
        if numStartPoints>0
            I = kmeansRounding(Y,numStartPoints,k,kmeansDistance);
            % compute content of each partition
            content = zeros(k,2);
            content(:,1) = histc(I,1:k);
            d = diag(D);
            for j = 1:k
                content(j,2) = sum(d(I==j));
            end
        end
        
    case 'graclus'

        I = graclusRounding(W,gamma,Pi,k);
        % compute content of each partition
        content = zeros(k,2);
        content(:,1) = histc(I,1:k);
        d = diag(D);
        for j = 1:k
            content(j,2) = sum(d(I==j));
        end
end

end

function [W,k,Y0,options,defaultopt] = parseInputs(varargin)
% parse input arguments

% check number of inputs
narginchk(1,4);

% get/check W
W = varargin{1};
if ~isreal(W) || ~isequal(ndims(W),2) || ~isequal(size(W,1),size(W,2))
    error([mfilename,':InvalidW'],...
        'W must be real square weight matrix.');
end

% get/check k
if nargin<2
    k = [];
else
    k = varargin{2};
end
if isempty(k)
    k = 2;
end
if ~isscalar(k) || ~isequal(k,round(k)) || k<1
    error([mfilename,':InvalidK'],...
        'k must be positive integer.');
end

% get/check Y0
if nargin<3
    Y0 = [];
else
    Y0 = varargin{3};
end
if ~isempty(Y0)
    if ~isreal(Y0) || ~isequal(ndims(Y0),2) || ~isequal(size(Y0),[size(W,1),k-1])
        error([mfilename,':InvalidY0'],...
            'Y0 must be real size(W,1) x (k-1) matrix.');
    end
end

% get options
if nargin<4
    options = [];
else
    options = varargin{4};
end
if isempty(options)
    options = ccbcutSet;
end

% set default options
defaultopt = struct('Display','final','Algorithm','irrq','Tau',1,...
    'BalanceType','normalized','RoundingType','kmeans',...
    'KmeansNumStartPoints',10,'KmeansDistance','sqeuclidean',...
    'MaxSOCIters',100,'MaxL1ProbIters',100,'TolBregSOCX',1e-3,...
    'TolBregSOCF',1e-3,'TolBregL1ProbX',1e-3,'TolBregL1ProbF',1e-3,...
    'BregLambda',1000,'R',100,'MinIRRQIters',5,'MaxIRRQIters',100,'TolIRRQX',1e-3,'TolIRRQF',1e-3,...
    'TolIRRQL1',1e-3,'TolIRRQL1Curv',1e-1,'Epsilon',1,'KRatio',0.5,'IRRQConstraint','none',...
    'IRRQLambda',1,'IRRQq',1,'IRRQAlpha',0.5,'EpsilonDenom',1,...
    'KRatioDenom',0.5); 

end

