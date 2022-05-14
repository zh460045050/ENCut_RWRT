function [ d, v ] = doRWspb( W, times, k )
%DONCUT 此处显示有关此函数的摘要
%   此处显示详细说明


%loop_weights = calSelfLoop(img);
%loop_idx = sub2ind([N, N], [1:N]', [1:N]');


new_W = W^times;
W(W~=0) = new_W(W~=0);

[N, ~] = size(W);

% new_W = W^times;
% W(W ~= 0) = new_W(W~=0);


D_inv = sparse(1:N,1:N, 1./sum(W));
W = D_inv * W;
I = sparse(1:N, 1:N, ones(N, 1));
L = I - W;

[d, v] = eigs(L + eps * speye(size(L)), k, 'sm');
end

