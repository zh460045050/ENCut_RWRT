function [ clustering ] = initialKMeans( img, k )
%INITIALNCUT 此处显示有关此函数的摘要
%   此处显示详细说明
[X, Y, Z] = size(img);
img_val = double(reshape(img, X*Y, Z));
clustering = kmeans(img_val,k);

end

