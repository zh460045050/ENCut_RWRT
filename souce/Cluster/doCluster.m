function [ labels ] = doCluster( img, k, type )
%DOCLUSTER 此处显示有关此函数的摘要
%   此处显示详细说明

if length(type) == 5 && sum(type == 'RWCut') == 5
    labels = doRWNCutCluster(img, 5, k);
elseif length(type) == 3 && sum(type == 'Cut') == 3
    labels = doCutCluster(img, k);
elseif length(type) == 4 && sum(type == 'NCut') == 4
    labels = doNCutCluster(img, k);
elseif length(type) == 5 && sum(type == 'WNCut') == 5
    labels = doWNCutCluster(img, k);
elseif length(type) == 6 && sum(type == 'WRWCut') == 6
    labels = doWRWCutCluster(img, 5, k);
elseif length(type) == 7 && sum(type == 'MERWCut') == 7
    
elseif length(type) == 6 && sum(type == 'CCBCut') == 6
    labels = doCCBCutCluster( img, k );
elseif length(type) == 8 && sum(type == 'CCBRWCut') == 8
    
end


end

