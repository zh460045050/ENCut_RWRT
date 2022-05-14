function [ E_RW ] = calRWEnergy( A, labels, i )
%CALRWENERGY 此处显示有关此函数的摘要
%   此处显示详细说明

%%%LRW%%%
    [X, Y] = size(labels);
    N = X * Y;
    S_t = double(labels == i);
    I = sparse(1:N,1:N,ones(N,1)); 
    lines = zeros(N,1);
    label_idx = find(S_t(:)==1);
    
    
    nb = calNB(labels, i);
    eg = [];
    for j = 1:length(nb)
        [ieg_x, ieg_y] = find(labels == nb(j));
        ieg = sub2ind([X, Y], ieg_x, ieg_y);
        eg = [eg;ieg];
    end
    
    W = sparse(1:N, 1:N, zeros(N, 1));
    W(eg, eg) = A(eg, eg);
    %A = W;
    
    D_inv = sparse(1:N,1:N,1./sum(W)); % get 
    Mk = size(label_idx,1);
    lines(label_idx,1) = 1/Mk;%average of probabilities,lines=1/Mk*blk
    clear label_idx;
    S = D_inv * W; 
    cmtime_cmp=(I-0.990*S) \ lines;
    %E_RW = cmtime_cmp(:,1) ./ Mk;
    E_RW = cmtime_cmp(:,1);
    
    
end

