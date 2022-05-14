function [ labels ] = fastENCut( img, times, k, numDownsamp )
%DOWNSAMPCLUSTER 此处显示有关此函数的摘要
%   此处显示详细说明

salmap = doSaliency(img);

img = im2double(img);
sigma = 50;
Bs = cell(numDownsamp,1);
[X, Y, Z] = size(img);
N = X * Y;
imgVals = reshape(img,N,Z); clear Lab;
son_idx = [1:N]';
[ edges ] = createEdges( son_idx, X, Y);
edges = [edges;[1:N; 1:N]'];
weights = makeweights_old(edges,imgVals,sigma);
W = adjacency(edges,weights,N); clear edges weights;
W = (W + W')/2;

W = (W - min(min(W))) / (max(max(W)) - min(min(W)));

slf_ind = sub2ind(size(W), [1:N]', [1:N]');
W(slf_ind) = 1 * mapminmax( (mapminmax(full(sum(W)), 0, 1)' .* reshape(salmap, [N, 1]))', 0, 1);

Wd = W;

imSizeD = [X, Y];

Bs = cell(numDownsamp,1);
B = sparse(1:N, 1:N, ones(N, 1));
for i = 1 : numDownsamp
    [Nd, Nd] = size(Wd);
    [r,c] = ind2sub(imSizeD, 1:size(Wd,1));
    do_keep = (mod(r, 2) == 0) & (mod(c, 2) == 0);
    B = Wd(:, do_keep);
    Wd = B' * B;
    Bs{i} = B;
    aa = floor(times / 2);
    if mod(aa, 2) == 0
        aa = aa + 1;
    end
    times = aa;
end

Wd = Wd^(times);

if numDownsamp ~= 0
    B = Bs{1};
end
for i = 2:numDownsamp
    B = B * Bs{i};
end
Wd = B * Wd * B';

W(W~=0) = Wd(W~=0);


W(slf_ind) = 0;

W = (W + W')/2;


D_inv = sparse(1:N,1:N, 1./sum(W));
I = sparse(1:N, 1:N, ones(1, N));
W = D_inv * W;
L = I - W;



[d, v] = eigs(L + 10^(-10) * speye(size(L)), k, 'sm');

v = sum(v);
d = d(:,2:end);
v = v(:,2:end);
d = d ./ sqrt(v);

idx = kmeans(d, k ,'Distance', 'sqeuclidean',  'MaxIter', 1000, 'Replicates', 10);

labels = reshape(idx, [X, Y]);

end