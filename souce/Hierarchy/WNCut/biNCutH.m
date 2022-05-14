function [ V, b, test_lab, ncut ] = biNCutH( X, preV, sigma, beta, lamda, max_iter, max_bfgs_iter )
%BINCUTH 此处显示有关此函数的摘要
%   此处显示详细说明

[n, d] = size(X); %获取样本数及维度
V = preV; %初始化映射矩阵v
cost_list = [];
iter = 1;
v_list = [];
b_list = [];
count = 0;
%开始最优化过程
while iter < max_iter
%给定V求最优化B
    normV = V / norm(V, 1); %归一化矩阵V
    nX = X * normV'; %映射到新空间
    [sortnX, ~] = sort(nX); %从小到大排序
    [NCut_k, K_cut, K_ass_l, K_ass_r] = mex_optimalb(sortnX, sigma); %求每个位置的NCut值
    %[NCut_k, K_cut, K_ass_l, K_ass_r] = calNCut_k(sortnX, sigma);
    [ncut, k] = min(NCut_k); %求最优b值
    %ncut
    b = (sortnX(k) + sortnX(k+1)) / 2;
    
%求错误率
    %cost = 1 - max(sum(test_lab == lab), sum(-test_lab == lab)) / length(lab);
    %cost_list(iter) = cost;
    
    %{
    if cost < beta
        break;
    end
    %}
    
%给定b求最优V,使用bfgs算法
    
    %newV = bfgs( V, X, max_bfgs_iter, lamda, sigma);
    [newV, b, ncut] = mex_bfgs(V, X, max_bfgs_iter, lamda, sigma);
    iter = iter + 1;
%     newV = newV ./ norm(newV, 1);
%     if sum(abs(newV - V)) < 1e-5
%         iter 
%         break;
%     end
%     
    V = newV;
    
    %[ dV ] = cal_dV( V, X, k, sigma , K_cut, K_ass_l, K_ass_r ); %求NCut对V的导数
    %V = V - 1e-5 * dV;
end

normV = V / norm(V, 1); %归一化矩阵V
nX = X * normV'; %映射到新空间
[sortnX, ~] = sort(nX); %从小到大排序
%[NCut_k, K_cut, K_ass_l, K_ass_r] = calNCut_k(sortnX, sigma); %求每个位置的NCut值
%[ncut, k] = min(NCut_k); %求最优b值
%ncut
[NCut_k] = mex_optimalb(sortnX, sigma); %求每个位置的NCut值
[ncut, k] = min(NCut_k); %求最优b值


b = (sortnX(k) + sortnX(k+1)) / 2;
    
%求聚类标签
normV = V / norm(V, 1);
test_lab = zeros([n,1]);
test_lab(X * normV' - b >0) = 1;
test_lab(test_lab ~= 1) = -1;

%画图
%figure();
%plot([1:length(cost_list)], cost_list);
end

