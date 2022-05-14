function [ labels ] = doRWNCutCluster( img, times, k, sigma )
%DONCUT 此处显示有关此函数的摘要
%   此处显示详细说明
if nargin == 3
    sigma = 50;

img = im2double(img);

[X, Y, Z] = size(img);
N = X * Y;
imgVals = reshape(img,N,Z); clear Lab;
son_idx = [1:N]';
edges = createEdges( son_idx, X, Y );
weights = makeweights_old(edges,imgVals,sigma);
W = adjacency(edges,weights,N); clear edges weights;

W = (W + W') / 2;


new_W = W^times;
W(W~=0) = new_W(W~=0);
D_inv = sparse(1:N,1:N, 1./sum(W));
W = D_inv * W;
I = sparse(1:N, 1:N, ones(1,N));
L = I - W; 
[d, v] = eigs(L + (10^-10) * speye(size(L)), k, 'sm');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Cluster
% for i = 2:k
%     d_i = d(:, i);
%     d_i = mapminmax(d_i', 0, 1)';
%     d_i = reshape(d_i, [X, Y]);
%     figure();imshow(d_i);
%     
% %     d_i = d(:, i);
% %     idx = d_i;
% %     idx(idx < 0 ) = 0;
% %     idx(idx ~= 0 ) = 1;
% %     idx = 1-idx;
% %     labels = reshape(idx, [X, Y]);
% %     result = show_result(img, labels);
% %     figure();imshow(result);
%         
% end
%  
% 
% d = d(:, 2:k);

%d = mapminmax(d', 0, 1)';


%idx = kmeans(d, k, 'Replicates', 10, 'MaxIter', 1000);
%d = d(:,2:end);

idx = kmeans(d, k ,'Distance', 'sqeuclidean', 'Replicates', 10, 'MaxIter', 1000);

labels = reshape(idx, [X, Y]);

end

