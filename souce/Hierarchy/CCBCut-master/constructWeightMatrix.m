function W = constructWeightMatrix(img,neighbors,beta)
% construct weight matrix for image

% author: Nathan D. Cahill
% email: nathan.cahill@rit.edu
% date: 21 February 2018

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
if nargin<3
    beta = 90;
end
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