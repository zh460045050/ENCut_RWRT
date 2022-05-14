function [ clustering] = initialNCut( W, k )
%INITIALNCUT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[N,~] = size(W);
D_inv = sparse(1:N,1:N, 1./sum(W));
W = D_inv * W;
I = sparse(1:N, 1:N, ones(N, 1));
L = I - W;
[d, v] = eigs(L + 10^(-10) * speye(size(L)), k, 'sm');
v = sum(v);
d = d(:,2:end);
v = v(:,2:end);
d = d ./ sqrt(v);
clustering = kmeans(d, k ,'Distance', 'sqeuclidean', 'Replicates', 10, 'MaxIter', 1000);
end

