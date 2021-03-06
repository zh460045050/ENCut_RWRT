function [ d, son_idx, ncut ] = biNCut( img, root_labels, lab)
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
if N < (X*Y) / 100.0
    labels = root_labels;
    ncut = 100;
    son_idx = -1;
    d = -1;
    normV = -1;
    return;
end


edges = createEdges( son_idx, X, Y );
weights = makeweights_old(edges,sonVals,sigma); %NCut
W = adjacency(edges,weights,N); clear edges weights;

D = sparse(1:N,1:N, sum(W)); % get 
D_inv = sparse(1:N,1:N, 1./sum(W));
%L = iD_inv * (D - W) * iD_inv;
L = D - W;
[d, v] = eigs(L, 2, 'sm');
v = sum(v);
[~, sec] = max(v);
d = d(:,sec);
labels = d;


max_lab = max(max(root_labels));
labels(labels < 0 ) = lab;
labels(labels ~= lab ) = max_lab + 1;

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

