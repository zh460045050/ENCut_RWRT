function [ normV, ncut ] = getV( X, preV )
%WCUT 此处显示有关此函数的摘要
%   此处显示详细说明


%NCutH求V
sigma_H = 1e-1;
beta = 1e-5;
lamda = 1e-1;
max_iter = 2;
max_bfgs_iter = 100;

[ normV, ~, ~, ncut ] = biNCutH( X, preV, sigma_H, beta, lamda, max_iter, max_bfgs_iter );

normV = normV / norm(normV ,1);
end

