function [ labels ] = dCut( img, k, numDownsamp )
%DOWNSAMPCLUSTER 此处显示有关此函数的摘要
%   此处显示详细说明
img = im2double(img);
sigma = 50;
Bs = cell(numDownsamp,1);
[X, Y, Z] = size(img);
N = X * Y;
imgVals = reshape(img,N,Z); clear Lab;
son_idx = [1:N]';
[ edges ] = createEdges( son_idx, X, Y);
weights = makeweights_old(edges,imgVals,sigma);
W = adjacency(edges,weights,N); clear edges weights;
W = (W + W')/2;

Wd = W;


imSizeD = [X, Y];

for i = 1:numDownsamp
  
  % Create a binary array of the pixels that will remain after decimating
  
  % every other row and column
  [r,c] = ind2sub(imSizeD, 1:size(Wd,1));
  do_keep = (mod(r, 2) == 0) & (mod(c, 2) == 0);
  
  % Downsample the affinity matrix
  Ws = Wd(:,do_keep)';
  
  % Downsample Y0 if provided on finest level
  
  % Normalize the downsampled affinity matrix
  d = (sum(Ws,1) + eps);
  B = bsxfun(@rdivide, Ws, d)';
  
  % "Square" the affinity matrix, while downsampling
  Wd = Ws*B;

  imSizeD = floor(imSizeD / 2);
  
  % Hold onto the normalized affinity matrix for bookkeeping
  Bs{i} = B;  

end


%Wd = (Wd + Wd')/2;

[Nd, Nd] = size(Wd);
Ds = cell(numDownsamp+1,1);

D = sparse(1:Nd,1:Nd, sum(Wd)); % get 
D_inv = sparse(1:Nd,1:Nd, 1./sum(Wd));
iD_inv=sqrt(D_inv);
L = D - Wd;
[Ds{numDownsamp+1}, v] = eigs(L, k, 'sm');



for i = numDownsamp:-1:1
  Ds{i} = Bs{i} * Ds{i+1};
end


d = Ds{1};

idx = kmeans(d, k ,'Distance', 'sqeuclidean', 'Replicates', 10, 'MaxIter', 1000);

labels = reshape(idx, [X, Y]);

end

