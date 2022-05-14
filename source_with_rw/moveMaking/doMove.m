function [ prob, unaries ] = KernelBound( A, K, current_clustering, energy_type)
%KERNELBOUND derives linear upper bound w.r.t binary indicator variables
%for NC (Normalized Cut) or AA (Average Association)
%Inputs:
%   A --- Affinity matrix (sparse or full) of N-by-N
%   K --- Number of desired clusters, default is 2
%   current_clustering --- N-by-1 vector of cluster indicator, value from 1
%   to K, default is a random generated vector
%   energy_type --- String 'NC' or 'AA', default is 'NC'
%Outputs:
%   unaries --- N-by-K matrix where unaries(n,k) is the unary (linear) cost
%   of assigning data n to cluster k


N = size(A,1); % number of data points

unaries = zeros(N, K, 'double');
% degree vector
d = sum(A,2);
%figure();imshow(reshape(mapminmax(d',0,1)', 321, 481));
prob = zeros(N, K);

for i=1:K
    % current binary indicators (N-by-1) for cluster i
    S_t = double(current_clustering == i);
    % compute gradient as unaries (see kernel bound proposition in the paper)
    if strcmp(energy_type, 'NC') % for normalized cut
        E_Ncut = S_t'* A * S_t / (d' * S_t)^2 * d - 2 * A * S_t / (d' * S_t);
        unaries(:,i) = E_Ncut;
    elseif strcmp(energy_type, 'NCRW')
        E_RW  = calRWEnergy( A, current_clustering, i );
        E_Ncut = S_t'* A * S_t / (d' * S_t)^2 * d - 2 * A * S_t / (d' * S_t);
        
        %E_RW = (E_RW - min(E_RW)) / (max(E_RW) - min(E_RW));  
        
        %unaries(:,i) = E_Ncut - 1e-2 * (1 - mapminmax(d', 0, 1)') .* E_RW;
        %unaries(:,i) = E_Ncut - 10 * (1 - mapminmax(d', 0, 1)' + 1e-5) .* E_RW;
        unaries(:,i) = E_Ncut - 1e-2 * (1 - mapminmax(d', 0, 1)' + 1e-5) .* E_RW;
        %unaries(:,i) = E_Ncut - E_RW;
        prob(:, i) = E_RW;
    elseif strcmp(energy_type, 'RW')
        E_RW  = calRWEnergy( A, current_clustering, i );
        unaries(:,i) = -E_RW;
        prob(:, i) = E_RW;
    end
end
end

