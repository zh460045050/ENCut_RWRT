function [ nb ] = calNB( labels, i )
%CALNB 此处显示有关此函数的摘要
%   此处显示详细说明

[X, Y] =size(labels);

[x, y] = find(labels == i);
idx = sub2ind([X, Y], x, y);

nb_idx = [];
for i = -1:1
    for j = -1:1
        inb_idx = idx + i * X + j;
        nb_idx = [nb_idx; inb_idx];
    end
end
nb_idx(nb_idx > X*Y) = [];
nb_idx(nb_idx <= 0) = [];
nb = unique(labels(nb_idx));
end
