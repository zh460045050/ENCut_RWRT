function [labels,p] = kmeansRounding(varargin)
% kmeansRounding: estimate labels from generalized eigenvectors
% usage: labels = kmeansRounding(V);
%    or: labels = kmeansRounding(V,numStartPoints);
%    or: labels = kmeansRounding(V,numStartPoints,k);
%    or: labels = kmeansRounding(V,numStartPoints,k,kmeansDistance);
%    or: [labels,p] = kmeansRounding(...);
%
% arguments:
%   V (nxm) - m generalized eigenvectors 
%   numStartPoints - number of times to run algorithm from random starting
%       points. Default = 10. K-means rounding will be run repeatedly,
%       and the mode for each point will be assigned as the class label
%   k (scalar) - number of clusters. Default k = m+1.
%   kmeansDistance (string) - distance measure used by kmeans function. Can
%       be one of 'sqeuclidean' (default), 'cityblock', 'cosine',
%       'correlation', or 'hamming'
%
%   labels (nx1) - class labels for each row of V
%   p (nx1) - mode frequencies
%

% author: Nathan D. Cahill
% email: nathan.cahill@rit.edu
% date: 21 February 2018

% parse input arguments
[V,numStartPoints,k,kmeansDistance] = parseInputs(varargin{:});

% initialize labels for multiple trials
labelTrials = zeros(size(V,1),numStartPoints);

% perform K-means clustering multiple times
for i = 1:numStartPoints
    labelTrials(:,i) = kmeans(V,k,'Distance',kmeansDistance);
end

% loop through each trial, permuting labels to best match initial trial
for i = 2:numStartPoints
    labelTrials(:,i) = bestPermutation(labelTrials(:,1),labelTrials(:,i));
end

% compute mode across trials
[labels,p] = mode(labelTrials,2);

% frequency of each mode
p = p./numStartPoints;

end

function [V,numStartPoints,k,kmeansDistance] = parseInputs(varargin)
% parse input arguments

% check number of inputs
narginchk(1,4)

% get/check V
V = varargin{1};
m = size(V,2);
if ~isequal(ndims(V),2)
    error([mfilename,':InvalidV'],...
        'V must be nxk array.');
end

% get/check numStartPoints
if nargin<2
    numStartPoints = [];
else
    numStartPoints = varargin{2};
end
if isempty(numStartPoints)
    numStartPoints = 10;
end
if ~isscalar(numStartPoints) || ~isequal(numStartPoints,round(numStartPoints)) || numStartPoints<1
    error([mfilename,':InvalidNumStartPoints'],...
        'numStartPoints must be positive integer.');
end

% get/check k
if nargin<3
    k = [];
else
    k = varargin{3};
end
if isempty(k)
    k = m+1;
end
if ~isscalar(k) || ~isequal(k,round(k)) || k<2
    error([mfilename,':InvalidK'],...
        'k must be integer greater than 1.');
end

% get/check kmeansDistance
if nargin<4
    kmeansDistance = '';
else
    kmeansDistance = varargin{4};
end
if isempty(kmeansDistance)
    kmeansDistance = 'sqeuclidean';
end
switch kmeansDistance
    case {'sqeuclidean','cityblock','cosine','correlation','hamming'}
    otherwise
        error([mfilename,':InvalidKmeansDistance'],...
            'kmeansDistance must be one of ''sqeuclidean'', ''cityblock'', ''cosine'', ''correlation'', or ''hamming''.');
end

end