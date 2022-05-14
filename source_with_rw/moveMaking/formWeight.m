function [ W ] = formWeight( img, sigma, type )
%CALENCUTWEIGHT 此处显示有关此函数的摘要
%   此处显示详细说明
if nargin < 3
    type = 'ENCut';
end
if nargin < 2
    sigma = 50;
end

times = 5;
img = im2double(img);
[X, Y, Z] = size(img);
N = X * Y;
son_idx = [1:N]';
imgVals = reshape(img,N,Z); clear Lab;
edges = createEdges( son_idx, X, Y );
weights = makeweights_old(edges,imgVals,sigma);
W = adjacency(edges,weights,N); clear edges weights;
W = (W + W') / 2;

if strcmp(type, 'ENCut')
    W = (W - min(min(W))) / (max(max(W)) - min(min(W)));
    W = (W + W') / 2;
    new_W = W^times;
    avg_W = new_W;
    avg_W(avg_W ~= 0) = 1;
    d_avg = sum(avg_W);
    W(W~=0) = new_W(W~=0);
    W = W ./ d_avg ;
    W = (W - min(min(W))) / (max(max(W)) - min(min(W)));
    W = (W + W') / 2;
    slf_ind = sub2ind(size(W), [1:N]', [1:N]');
    W(slf_ind) = 0;
end

W = (W - min(min(W))) / (max(max(W)) - min(min(W)));

end

