function [ d, son_idx, ncut ] = biSNCut( img, times, root_labels, lab, salmap)
%WNCUT 此处显示有关此函数的摘要
%   此处显示详细说明
%   chosen = 1 使用WNCut
%   chosen = 0 使用NCut

%初始化参数
[X, Y, Z] = size(img);
r = img(:,:,1);
g = img(:,:,2);
b = img(:,:,3);
sigma = 50;



son_idx = sub2ind([X, Y], find(root_labels == lab));
son_img = [r(son_idx), g(son_idx), b(son_idx)];
N = length(son_idx);
sonVals = reshape(son_img, N, Z);
% if N < (X*Y) / 100.0
%     labels = root_labels;
%     ncut = 100;
%     son_idx = -1;
%     d = -1;
%     normV = -1;
%     return;
% end

cur_salmap = salmap(son_idx);

edges = createEdges( son_idx, X, Y );
edges = [edges;[1:N; 1:N]'];
weights = makeweights_old(edges,sonVals,sigma); %NCut
W = adjacency(edges,weights,N); clear edges weights;

W = (W + W') / 2;


W = (W - min(min(W))) / (max(max(W)) - min(min(W)));

slf_ind = sub2ind(size(W), [1:N]', [1:N]');
W(slf_ind) = 1 * mapminmax( (mapminmax(full(sum(W)), 0, 1)' .* reshape(cur_salmap, [N, 1]))', 0, 1);


new_W = W^times;
avg_W = new_W;
avg_W(avg_W ~= 0) = 1;
d_avg = sum(avg_W);
W(W~=0) = new_W(W~=0);
W = W ./ d_avg;

%W(sub2ind([N, N], [1:N]', [1:N]')) = 0;

W(slf_ind) = 0;

W = (W + W') / 2;


W = (W - min(min(W))) / (max(max(W)) - min(min(W)));

slf_ind = sub2ind(size(W), [1:N]', [1:N]');
W(slf_ind) = 0;

org_W = W;

D_inv = sparse(1:N,1:N, 1./sum(W));
I = sparse(1:N, 1:N, ones(1,N));
W = D_inv * W;
L = I - W;
%L = D - W;
[d, v] = eigs(L + 10^(-10) * speye(size(L)), 2, 'sm');
v = sum(v);
[~, sec] = max(v);
d = d(:,sec);
labels = d;

W = org_W;

max_lab = max(max(root_labels));


labels(labels < 0 ) = lab;
labels(labels ~= lab ) = max_lab + 1;

% idx = kmeans(d, 2, 'Distance', 'sqeuclidean',  'MaxIter', 10000, 'Replicates', 20);
% labels(idx == 1) = lab;
% labels(idx ~= lab) = max_lab + 1;


root_labels(son_idx) = labels;
labels = root_labels;



turn_son_idx = zeros(N, 1) - 1;
turn_son_idx(son_idx([1:N]')) = [1:N]';


lab_l = lab;
son_idx_l = sub2ind([X, Y], find(labels == lab_l));
edges_l = createEdges( son_idx_l, X, Y );
idx_l = sub2ind(size(W), turn_son_idx(son_idx_l(edges_l(:,1))), turn_son_idx(son_idx_l(edges_l(:,2))));
lab_r = max(max(root_labels));
son_idx_r = sub2ind([X, Y], find(labels == lab_r));
edges_r = createEdges( son_idx_r, X, Y );
idx_r = sub2ind(size(W), turn_son_idx(son_idx_r(edges_r(:,1))), turn_son_idx(son_idx_r(edges_r(:,2))));
d_l = sum(sum(W)) - sum(W(idx_l)); 
d_r = sum(sum(W)) - sum(W(idx_r));
cut = sum(sum(W)) - sum(W(idx_l)) - sum(W(idx_r));
ncut = cut / d_l + cut / d_r;

end

