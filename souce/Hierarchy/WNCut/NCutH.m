function [ final_lab, line_list ] = NCutH( X, num_cluster, sigma, beta, max_iter, lamda, max_bfgs_iter )
%%%%%%%%%%%%%%%%%%%%%%%%
%   X:��������
%   num_cluster:�����
%   sigma:��˹�������
%   beta:�����ʧ������������
%   lamba:bfgs�㷨����
%   max_iter:��������������
%   max_bfgs_iter:bfgs����������
%%%%%%%%%%%%%%%%%%%%%%%%
%       ������������
%%%%%%%%%%%%%%%%%%%%%%%%
%�������
if nargin < 7
    max_bfgs_iter = 1000;
end
if nargin < 6
    lamda = 1e-3;
end
if nargin < 5
    max_iter = 100;
end
if nargin < 4
    beta = 1e-5;
end
if nargin < 3
    sigma = 1;
end

%��ʼ��
life_node = 1;
count_lab = 0;
[n, d] = size(X);
final_lab = zeros(n, 1);
line_list = {};
count_line = 0;
queue = {};


%���е�һ�ξ���
[ V, b, test_lab, ncut ] = biNCutH( X, sigma, beta, lamda, max_iter, max_bfgs_iter);

%�����������
queue{1, 1} = -ncut; 
queue{1, 2} = X;
queue{1, 3} = V;
queue{1, 4} = b;
queue{1, 5} = test_lab;
queue{1, 6} = 1:n;
queue{1, 6} = queue{1, 6}';
count_que = 1;

%�������о���ֱ���������ո���
while life_node < num_cluster
    %ȡ��NCutֵ��С��Ҷ�ӽڵ㣬���Ӷ�����ɾ��
    cur_X = queue{count_que, 2};
    cur_V = queue{count_que, 3};
    cur_b = queue{count_que, 4};
    cur_test_lab = queue{count_que, 5};
    cur_idx = queue{count_que, 6};
    
    count_line = count_line + 1;
    line_list{count_line}.V = cur_V;
    line_list{count_line}.b = cur_b;
    
    queue(count_que, :) = [];
    count_que = count_que - 1;
    
    %�Ըýڵ��и����������б��
    count_lab = count_lab + 1;
    final_lab(cur_idx(cur_test_lab == 1)) = count_lab;
    life_node = life_node + 1;
    
    if life_node == num_cluster
        break;
    end
    
    
    %������������и��NCut����
    l_idx = cur_idx(sub2ind(size(cur_X), find(cur_test_lab == 1)));
    l_X = X(l_idx, :);
    if length(l_idx) > 1
        [ V_l, b_l, test_lab_l, ncut_l ] = biNCutH( l_X, sigma, beta, lamda, max_iter, max_bfgs_iter);
        count_que = count_que + 1;
        queue{count_que, 1} = -ncut_l; 
        queue{count_que, 2} = l_X;
        queue{count_que, 3} = V_l;
        queue{count_que, 4} = b_l;
        queue{count_que, 5} = test_lab_l;
        queue{count_que, 6} = l_idx;
    end
    %�������Ҷ����и��NCut����
    r_idx = cur_idx(sub2ind(size(cur_X), find(cur_test_lab == -1)));
    r_X = X(r_idx, :);
    if length(r_idx) > 1
        [V_r, b_r, test_lab_r, ncut_r ] = biNCutH( r_X, sigma, beta, lamda, max_iter, max_bfgs_iter);
        count_que = count_que + 1;
        queue{count_que, 1} = -ncut_r; 
        queue{count_que, 2} = r_X;
        queue{count_que, 3} = V_r;
        queue{count_que, 4} = b_r;
        queue{count_que, 5} = test_lab_r;
        queue{count_que, 6} = r_idx;
    end

    %���и���Ncutֵ��������
    queue = sortrows(queue, 1);
end

%ת��Ϊ0�ı�ǩ
final_lab = final_lab + 1;

end

