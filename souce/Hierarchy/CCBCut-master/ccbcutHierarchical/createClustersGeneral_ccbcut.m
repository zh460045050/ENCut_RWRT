function [allClusters, cut, cutPart1, cutPart2, threshold] =  createClustersGeneral_ccbcut(vmin,W,opts,threshold_type,deg)
% Transforms an eigenvector into a cluster indicator function by thresholding.
% 
% Usage: [allClusters, cut, cutPart1, cutPart2, threshold] 
%			= createClusters_ccbcut(vmin,W,opts,threshold_type);
%
% Input:
%   vmin: The eigenvector.
%   W: Weight matrix.
%   opts: CCBCut options structure.
%   threshold_type: 0: zero, 1: median, 2: mean, -1: best
%
%   (optional) deg: Degrees of vertices as column vector. Default is 
%   sum(W,2) in normalized case and ones(size(W,1),1) in unnormalized
%   case. Will be ignored if normalized=false.
%
% Output:
%   allClusters: Obtained clustering after thresholding.
%   cut: Value of the CCBCut.
%   cutpart1,cutpart2: The two components of Ratio/Normalized Cut and 
%   Ratio/Normalized Cheeger Cut.
%   threshold: The threshold used to obtain the partitioning.
%
% author: Nathan D. Cahill, based on modifications of code from:
%
% (C)2010-11 Thomas Buehler and Matthias Hein
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
	
% Default values for deg
if (nargin<6)
    if isequal(ccbcutGet(opts,'BalanceType'),'normalized')
        deg=sum(W,2);
    else
        deg=ones(size(W,1),1);
    end
end
%Make deg a row vector;
if (size(deg,1)>1)
    deg=deg';
end

assert(isempty(find(diag(W)~=0,1)),'Graph contains self loops. W has to have zero diagonal.');

if threshold_type>=0
    threshold= determineThreshold(threshold_type,vmin);
    %allClusters= computeClusterIndicatorFunction(vmin,threshold);
    allClusters= (vmin>threshold);
    [cutPart1,cutPart2] = computeCutValue_ccbcut(allClusters,W,opts,deg); %cutPart1: vmin<threshold, cutPart2: vmin>threshold
    cut=cutPart1+cutPart2;
else
    
    [vmin_sorted, index]=sort(vmin);
    W_sorted=W(index,index);
    
    % sum of all degrees in the cluster minus weights within cluster
    deg2=sum(W_sorted); % this has to be the degree also in unnormalized variant
    tempcuts_threshold=cumsum(deg2) - 2*cumsum(full(sum(triu(W_sorted))));
    
    % divide by volume/size
    tau = ccbcutGet(opts,'Tau');
    if isequal(ccbcutGet(opts,'BalanceType'),'normalized')
        volumes_threshold=cumsum(deg(index));
        
        cutparts1_threshold=tempcuts_threshold(1:end-1)./((2*volumes_threshold(1:end-1)).^(tau/2));
        cutparts1_threshold(isnan(cutparts1_threshold))=0;
        cutparts2_threshold=tempcuts_threshold(1:end-1)./((2*(volumes_threshold(end)-volumes_threshold(1:end-1))).^(tau/2));
        cutparts2_threshold(isnan(cutparts2_threshold))=0;
    else
        sizes_threshold=cumsum(ones(1,size(vmin,1)-1));
        cutparts1_threshold=tempcuts_threshold(1:end-1)./((2*sizes_threshold).^(tau/2));
        cutparts2_threshold=tempcuts_threshold(1:end-1)./((2*(size(vmin,1)-sizes_threshold)).^(tau/2));
    end
    
    % calculate cuts/cheegers
    cuts_threshold=cutparts1_threshold+cutparts2_threshold;    
    
    % don't threshold within regions of same value
    
    [~,indexU]=unique(vmin_sorted(1:end-1));% unique gives index of last occurence
    %[~,indexU]=unique(vmin_sorted);% unique gives index of last occurence
    
    % find best cut/cheeger
    [cut,threshold_index]=min(cuts_threshold(indexU));
    
    % update
    cutPart1=cutparts1_threshold(indexU(threshold_index));
    cutPart2=cutparts2_threshold(indexU(threshold_index));
    
    threshold=vmin_sorted(indexU(threshold_index));
    allClusters=vmin>threshold;
    
end

end

function threshold = determineThreshold(threshold_type,u)
% Select the treshold for the cluster indicator function

assert(threshold_type==0 || threshold_type==1 || threshold_type==2);

switch threshold_type
    case 0
        threshold=0;
    case 1
        threshold = median(u);
    case 2
        threshold = mean(u);
end

end