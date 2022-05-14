function [ labels ] = doCutCluster( img, k, sigma )
%DONCUT 此处显示有关此函数的摘要
%   此处显示详细说明
if nargin == 2
    sigma = 50;

img = im2double(img);

[X, Y, Z] = size(img);
N = X * Y;
imgVals = reshape(img,N,Z); clear Lab;
son_idx = [1:N]';
edges = createEdges( son_idx, X, Y );
weights = makeweights_old(edges,imgVals,sigma);
W = adjacency(edges,weights,N); clear edges weights;

D = sparse(1:N,1:N, sum(W)); % get 



L = D - W;
[d, v] = eigs(L + (10^-10) * speye(size(L)), k, 'sm');

idx = kmeans(d, k);

labels = reshape(idx, [X, Y]);

end

