function [ labels,d  ] = doENCutSaliency( img, times, k, sigma )
%DONCUT 此处显示有关此函数的摘要
%   此处显示详细说明
if nargin == 3
    sigma = 50;

salmap = doSaliency(img);
img = im2double(img);

[X, Y, Z] = size(img);
N = X * Y;
son_idx = [1:N]';
imgVals = reshape(img,N,Z); clear Lab;
edges = createEdges( son_idx, X, Y );
edges = [edges;[1:N; 1:N]'];
weights = makeweights_old(edges,imgVals,sigma);
W = adjacency(edges,weights,N); clear edges weights;

W = (W + W') / 2;


% rr = reshape(mapminmax(sum(W), 0, 1), [X, Y]);
% new_d = mapminmax(sum(W), 0, 1) + 2 * (mapminmax(sum(W), 0, 1) .* reshape(salmap, [1, N]));
% new_d = mapminmax(new_d, 0, 1)';
% rr_new = reshape(new_d, [X, Y]);
% % 
% figure();imshow(rr);
% figure();imshow(salmap);
% figure();imshow(rr_new);

W = (W - min(min(W))) / (max(max(W)) - min(min(W)));

slf_ind = sub2ind(size(W), [1:N]', [1:N]');
W(slf_ind) = 1 * mapminmax( (mapminmax(full(sum(W)), 0, 1)' .* reshape(salmap, [N, 1]))', 0, 1);
%W(slf_ind) = (mapminmax(full(sum(W)), 0, 1)' .* reshape(salmap, [N, 1]));
%W(slf_ind) = new_d;
%W(slf_ind) = 1 + mapminmax(full(sum(W)), 0, 1)' + 1 * reshape(salmap, [N, 1]);


%W = (W - min(min(W))) / (max(max(W)) - min(min(W)));

W = (W - min(min(W))) / (max(max(W)) - min(min(W)));

W = (W + W') / 2;


new_W = W^times;
avg_W = new_W;
avg_W(avg_W ~= 0) = 1;
d_avg = sum(avg_W);
W(W~=0) = new_W(W~=0);
W = W ./ d_avg ;

%W(sub2ind([N, N], [1:N]', [1:N]')) = 0;


W = (W - min(min(W))) / (max(max(W)) - min(min(W)));

%W = W + org_W;

W = (W + W') / 2;


% d = new_W(sub2ind(size(W), [1:N]', [1:N]'));
% d = mapminmax(d', 0, 1)';
% rr = reshape(d, [X, Y]);
% figure();imshow(rr);



slf_ind = sub2ind(size(W), [1:N]', [1:N]');
W(slf_ind) = 0;


W = (W - min(min(W))) / (max(max(W)) - min(min(W)));

D = sparse(1:N,1:N, sum(W)); % get 
D_inv = sparse(1:N,1:N, 1./sum(W));

W = D_inv * W;
I = sparse(1:N, 1:N, ones(N, 1));
L = I - W;
[d, v] = eigs(L + 10^(-10) * speye(size(L)), k, 'sm');
v = sum(v);
d = d(:,2:end);
v = v(:,2:end);
d = d ./ sqrt(v);
% 
% for i = 2:k
%     d_i = d(:, i);
%     d_i = mapminmax(d_i', 0, 1)';
%     d_i = reshape(d_i, [X, Y]);
%     figure();imshow(d_i);
% end

idx = kmeans(d, k ,'Distance', 'sqeuclidean', 'Replicates', 10, 'MaxIter', 1000);
%[~, idx] = vl_kmeans(d', k);

labels = reshape(idx, [X, Y]);


end

