function [ edges ] = createEdges( son_idx, X, Y, r )
%CREATEEDGES 此处显示有关此函数的摘要
%   此处显示详细说明
if nargin < 4
    r = 10;
end

root_N = X * Y;

N = length(son_idx);
nb = [];
for i = -r:r
    for j = -r:r
        if i == 0 && j == 0
            continue;
        end
        if i^2 + j^2 <= r^2
            nb = [nb;i, j];
        end
    end
end

edges = [];
for i = 1:length(nb)
    edges_l = [son_idx, son_idx + X * nb(i,1) + nb(i, 2)];
    edges = [edges;edges_l];
end
edges = [edges; [edges(:,2), edges(:,1)]];
excluded=find((edges(:,1)>root_N)|(edges(:,1)<1)|(edges(:,2)>root_N)| ...
    (edges(:,2)<1));
if ~isempty(excluded)
    edges(excluded, :)=[]; 
end

turn_son_idx = zeros(root_N, 1) - 1;
turn_son_idx(son_idx([1:N]')) = [1:N]';

edges = [turn_son_idx(edges(:,1)), turn_son_idx(edges(:,2))];

excluded=find((edges(:,1) == -1) |(edges(:,2) == -1));
if ~isempty(excluded)
    edges(excluded,:)=[]; 
end

edges = unique(edges, 'rows');

end

