function [ d, son_idx, ncut ] = biCCBCut( img, root_labels, lab)
%BICCBCUT 此处显示有关此函数的摘要
%   此处显示详细说明
optsCCNCut = ccbcutSet('Display','iter',...
    'Algorithm','irrq',...
    'Tau',1,...
    'BalanceType','normalized',...
    'RoundingType','kmeans',...
    'KmeansNumStartPoints',1,...
    'MinIRRQIters',5,...
    'MaxIRRQIters',50,...
    'TolIRRQX',eps,...
    'TolIRRQF',eps,...
    'TolIRRQL1',1e-5,...
    'TolIRRQL1Curv',0.08,...
    'Epsilon',1,...
    'KRatio',0.0001,...
    'IRRQConstraint','hard');

[X, Y, Z] = size(img);
r = img(:,:,1);
g = img(:,:,2);
b = img(:,:,3);
sigma = 50;


imBw = zeros(size(root_labels));
imBw(root_labels == lab) = 1;
imLabel = bwlabel(imBw);                %对各连通域进行标记
stats = regionprops(imLabel,'Area');    %求各连通域的大小
area = cat(1,stats.Area);
[~, a_idx] = max(area);

son_idx = sub2ind([X, Y], find(imLabel == a_idx));
son_img = [r(son_idx), g(son_idx), b(son_idx)];
N = length(son_idx);
sonVals = reshape(son_img, N, Z);
if N < (X*Y) / 100.0
    ncut = 100;
    son_idx = -1;
    d = -1;
    return;
end


edges = createEdges( son_idx, X, Y );
weights = makeweights_old(edges,sonVals,sigma); %NCut
W = adjacency(edges,weights,N); clear edges weights;



[~, d] = myccbcut(W,2,[],optsCCNCut);
labels = d;


max_lab = max(max(root_labels));


labels(labels < 0 ) = lab;
labels(labels ~= lab ) = max_lab + 1;

root_labels(son_idx) = labels;
labels = root_labels;



turn_son_idx = zeros(N, 1) - 1;
turn_son_idx(son_idx([1:N]')) = [1:N]';


temp_labels = root_labels;
temp_labels(imLabel ~= a_idx) = -1;
temp_labels(imLabel == a_idx) = root_labels(imLabel == a_idx);

lab_l = lab;
son_idx_l = sub2ind([X, Y], find(temp_labels == lab_l));
edges_l = createEdges( son_idx_l, X, Y );
idx_l = sub2ind(size(W), turn_son_idx(son_idx_l(edges_l(:,1))), turn_son_idx(son_idx_l(edges_l(:,2))));
lab_r = max(max(root_labels));
son_idx_r = sub2ind([X, Y], find(temp_labels == lab_r));

edges_r = createEdges( son_idx_r, X, Y );
idx_r = sub2ind(size(W), turn_son_idx(son_idx_r(edges_r(:,1))), turn_son_idx(son_idx_r(edges_r(:,2))));
d_l = sum(sum(W)) - sum(W(idx_l)); 
d_r = sum(sum(W)) - sum(W(idx_r));
cut = sum(sum(W)) - sum(W(idx_l)) - sum(W(idx_r));
ncut = cut / d_l + cut / d_r;

%ncut = ncut / (sum(sum(W)) / N);


end

