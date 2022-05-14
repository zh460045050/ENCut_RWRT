function [ d, v ] = doRWspb( W, times, k, salmap )
%DONCUT 此处显示有关此函数的摘要
%   此处显示详细说明


[N, ~] = size(W);
slf_ind = sub2ind(size(W), [1:N]', [1:N]');
W(slf_ind) = 1 * mapminmax( (mapminmax(full(sum(W)), 0, 1)' .* reshape(salmap, [N, 1]))', 0, 1);

new_W = W^times;
W(W~=0) = new_W(W~=0);

W(slf_ind) = 0;

[N, ~] = size(W);


D_inv = sparse(1:N,1:N, 1./sum(W));
W = D_inv * W;
I = sparse(1:N, 1:N, ones(N, 1));
L = I - W;

[d, v] = eigs(L + eps * speye(size(L)), k, 'sm');
end

