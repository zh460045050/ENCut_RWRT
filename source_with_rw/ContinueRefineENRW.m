function [ labels ] = ENRW( img, k, clustering, maxiter, sigma )
%ENRW 此处显示有关此函数的摘要
%   此处显示详细说明

if nargin < 5
    sigma = 50;
end
if nargin < 4
    maxiter = 5;
end

times = 5;
img = im2double(img);
[X, Y, Z] = size(img);
N = X * Y;
son_idx = [1:N]';
imgVals = reshape(img,N,Z); clear Lab;
edges = createEdges( son_idx, X, Y );
weights = makeweights_old(edges,imgVals,sigma);
W = adjacency(edges,weights,N); clear edges weights;
W = (W + W') / 2;
W = (W - min(min(W))) / (max(max(W)) - min(min(W)));
A_RW = W;

new_W = W^times;
avg_W = new_W;
avg_W(avg_W ~= 0) = 1;
d_avg = sum(avg_W);
W(W~=0) = new_W(W~=0);
W = W ./ d_avg ;
W = (W - min(min(W))) / (max(max(W)) - min(min(W)));
W = (W + W') / 2;
slf_ind = sub2ind(size(W), [1:N]', [1:N]');
W(slf_ind) = 0;
W = (W - min(min(W))) / (max(max(W)) - min(min(W)));
A_NC = W;


itrnum = 1;

while 1
    
    itrnum
    
    
    unaries = zeros(N, k, 'double');
    %A_NC = A_RW;
    d_RW = sum(A_RW,2);
    d_NC = sum(A_NC,2);
    for i=1:k
        % current binary indicators (N-by-1) for cluster i
        S_t = double(clustering == i);
        E_RW  = calRWEnergy( A_RW, clustering, i );
        E_Ncut = S_t'* A_NC * S_t / (d_NC' * S_t)^2 * d_NC - 2 * A_NC * S_t / (d_NC' * S_t);
        %unaries(:,i) = E_Ncut -  1e-1 * (1-S_t) .* E_RW;%
        %unaries(:,i) = E_Ncut -  1 * (1 - mapminmax(d_RW', 0, 1)' + 1e-5) .* E_RW;
        unaries(:,i) = E_Ncut - 100 .* E_RW;
        
        %%%%%Debug%%%%
%         RR_NCut = E_Ncut; RR_NCut = mapminmax(RR_NCut', 0 ,1); RR_NCut = reshape(RR_NCut, 161, 241);
%         RR_RW = E_RW; RR_RW = mapminmax(-RR_RW', 0 ,1); RR_RW = reshape(RR_RW, 161, 241);
%         RR_Sum = E_RW; RR_Sum = mapminmax(unaries(:,i)', 0 ,1); RR_Sum = reshape(RR_Sum, 161, 241);
%         source_lab = reshape(S_t, 161, 241);
%         figure();imshow(source_lab);
%         figure();imshow(RR_NCut);
%         figure();imshow(RR_RW);
%         figure();imshow(RR_Sum);
        %%%%%Debug%%%%
        
    end
    
    [~, new_clustering] = min(unaries, [], 2); 
    
    %showFigs(img, reshape(new_clustering, [X, Y]));
    % check convergence
    if 0 == sum(new_clustering~=clustering)
        disp('converged!');
        break;
    else
        clustering = new_clustering;
    end
    itrnum = itrnum + 1;
    if itrnum > maxiter
        disp('max iteration');
        break;
    end
end

labels = reshape(clustering, X, Y);
showFigs(img, labels);
end

