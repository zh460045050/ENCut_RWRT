function [ normV, ncut ] = getV( X, preV )
%WCUT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��


%NCutH��V
sigma_H = 1e-1;
beta = 1e-5;
lamda = 1e-1;
max_iter = 2;
max_bfgs_iter = 100;

[ normV, ~, ~, ncut ] = biNCutH( X, preV, sigma_H, beta, lamda, max_iter, max_bfgs_iter );

normV = normV / norm(normV ,1);
end

