function [ labels ] = doWRWCutCluster(  img, times, k, sigma  )
%DOWRWCUTCLUSTER 此处显示有关此函数的摘要
%   此处显示详细说明
if nargin == 3
    sigma = 50;

img = im2double(img);

r = img(:,:,1);
g = img(:,:,2);
b = img(:,:,3);


[X, Y, Z] = size(img);
N = X * Y;
imgVals = reshape(img,N,Z); clear Lab;
son_idx = [1:N]';
son_img = [r(son_idx), g(son_idx), b(son_idx)];
sonVals = reshape(son_img, N, Z);
preV = ones(1, Z);

[normV, ~] = getV(sonVals, preV);


edges = createEdges( son_idx, X, Y );
weights = makeweights_wncut(edges, imgVals, sigma, normV);
W = adjacency(edges,weights,N); clear edges weights;

W = (W + W') / 2;
new_W = W^times;
W(W ~= 0) = new_W(W ~= 0);
D_inv = sparse(1:N,1:N, 1./sum(W));
W = D_inv * W;
I = sparse(1:N, 1:N, ones(1,N));
L = I - W; 
[d, v] = eigs(L + (10^-10) * speye(size(L)), k, 'sm');

%d = d(:,2:end);

idx = kmeans(d, k ,'Distance', 'sqeuclidean', 'Replicates', 10, 'MaxIter', 1000);

labels = reshape(idx, [X, Y]);

end

