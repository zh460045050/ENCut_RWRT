function [ V, b, test_lab, ncut ] = biNCutH( X, preV, sigma, beta, lamda, max_iter, max_bfgs_iter )
%BINCUTH �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

[n, d] = size(X); %��ȡ��������ά��
V = preV; %��ʼ��ӳ�����v
cost_list = [];
iter = 1;
v_list = [];
b_list = [];
count = 0;
%��ʼ���Ż�����
while iter < max_iter
%����V�����Ż�B
    normV = V / norm(V, 1); %��һ������V
    nX = X * normV'; %ӳ�䵽�¿ռ�
    [sortnX, ~] = sort(nX); %��С��������
    [NCut_k, K_cut, K_ass_l, K_ass_r] = mex_optimalb(sortnX, sigma); %��ÿ��λ�õ�NCutֵ
    %[NCut_k, K_cut, K_ass_l, K_ass_r] = calNCut_k(sortnX, sigma);
    [ncut, k] = min(NCut_k); %������bֵ
    %ncut
    b = (sortnX(k) + sortnX(k+1)) / 2;
    
%�������
    %cost = 1 - max(sum(test_lab == lab), sum(-test_lab == lab)) / length(lab);
    %cost_list(iter) = cost;
    
    %{
    if cost < beta
        break;
    end
    %}
    
%����b������V,ʹ��bfgs�㷨
    
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
    
    %[ dV ] = cal_dV( V, X, k, sigma , K_cut, K_ass_l, K_ass_r ); %��NCut��V�ĵ���
    %V = V - 1e-5 * dV;
end

normV = V / norm(V, 1); %��һ������V
nX = X * normV'; %ӳ�䵽�¿ռ�
[sortnX, ~] = sort(nX); %��С��������
%[NCut_k, K_cut, K_ass_l, K_ass_r] = calNCut_k(sortnX, sigma); %��ÿ��λ�õ�NCutֵ
%[ncut, k] = min(NCut_k); %������bֵ
%ncut
[NCut_k] = mex_optimalb(sortnX, sigma); %��ÿ��λ�õ�NCutֵ
[ncut, k] = min(NCut_k); %������bֵ


b = (sortnX(k) + sortnX(k+1)) / 2;
    
%������ǩ
normV = V / norm(V, 1);
test_lab = zeros([n,1]);
test_lab(X * normV' - b >0) = 1;
test_lab(test_lab ~= 1) = -1;

%��ͼ
%figure();
%plot([1:length(cost_list)], cost_list);
end

