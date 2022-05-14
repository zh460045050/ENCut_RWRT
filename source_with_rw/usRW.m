function [ labels ] = usRW( img, k, maxiter, sigma )
%ENRW 此处显示有关此函数的摘要
%   此处显示详细说明

if nargin < 4
    sigma = 50;
end
if nargin < 3
    maxiter = 20;
end

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

D_inv = sparse(1:N,1:N, 1./sum(W));
W = D_inv * W;
I = sparse(1:N, 1:N, ones(N, 1));
L = I - W;
[d, v] = eigs(L + 10^(-10) * speye(size(L)), k, 'sm');
v = sum(v);
d = d(:,2:end);
v = v(:,2:end);
d = d ./ sqrt(v);
clustering = kmeans(d, k ,'Distance', 'sqeuclidean', 'Replicates', 10, 'MaxIter', 1000);
init_lab = reshape(clustering, X, Y);
showFigs(img, init_lab);

itrnum = 1;

while 1
    itrnum
    unaries = zeros(N, k, 'double');
    for i=1:k
        % current binary indicators (N-by-1) for cluster i
        
        E_RW  = calRWEnergy( A_RW, clustering, i );
        unaries(:,i) = -E_RW;
        
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

