function [ edges ] = createEdges( son_idx, X, Y )
%CREATEEDGES 此处显示有关此函数的摘要
%   此处显示详细说明

root_N = X * Y;

N = length(son_idx);

edges_l = [son_idx, son_idx - 1];
edges_r = [son_idx, son_idx + 1];
edges_u = [son_idx, son_idx + X];
edges_d = [son_idx, son_idx - X];
edges_lu = [son_idx, son_idx - X - 1];
edges_ru = [son_idx, son_idx - X + 1];
edges_ld = [son_idx, son_idx + X - 1];
edges_rd = [son_idx, son_idx + X + 1];

edges = [edges_l;edges_r;edges_u;edges_d;edges_lu;edges_ru;edges_ld;edges_rd];
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

