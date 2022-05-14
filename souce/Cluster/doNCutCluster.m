function [ labels ] = doNCutCluster( img, k, sigma )
%DONCUT 此处显示有关此函数的摘要
%   此处显示详细说明
if nargin == 2
    sigma = 50;

img = im2double(img);

[X, Y, Z] = size(img);
N = X * Y;
son_idx = [1:N]';
imgVals = reshape(img,N,Z); clear Lab;
edges = createEdges( son_idx, X, Y );
weights = makeweights_old(edges,imgVals,sigma);
W = adjacency(edges,weights,N); clear edges weights;

W = (W + W') / 2;



D = sparse(1:N,1:N, sum(W)); % get 
D_inv = sparse(1:N,1:N, 1./sum(W));
iD_inv=sqrt(D_inv);
W = D_inv * W;
I = sparse(1:N, 1:N, ones(N, 1));
L = I - W;
[d, v] = eigs(L, k, 'sm');
% 
% for i = 2:k
%     d_i = d(:, i);
%     d_i = mapminmax(d_i', 0, 1)';
%     d_i = reshape(d_i, [X, Y]);
%     figure();imshow(d_i);
% end

idx = kmeans(d, k ,'Distance', 'sqeuclidean', 'Replicates', 10, 'MaxIter', 1000);

labels = reshape(idx, [X, Y]);


end

