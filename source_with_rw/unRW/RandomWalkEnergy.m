function [ rw_e ] = RandomWalkEnergy( A, clustering )
%NormalizedCutEnergy computes normalized cut given affinity A and labeling
%Inputs:
%   A --- Affinity matrix (sparse or full) of N-by-N
%   clustering --- N-by-1 vector of cluster indicator, start from 1
%Outputs:
%   ncut_e --- Normalized Cut energy

K = max(clustering); % number of clusters
d = sum(A,2); % degree vector
rw_e = 0;
for k=1:K
    cur_E  = calRWEnergy( A, clustering, k );
    rw_e = rw_e + sum(cur_E(clustering == k));
end
% ncut_e = K - nassoc_e

end