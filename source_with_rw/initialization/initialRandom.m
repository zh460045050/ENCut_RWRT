function [ clustering ] = initialRandom( img, k )
%INITIALNCUT 此处显示有关此函数的摘要
%   此处显示详细说明
[X, Y, Z] = size(img);
N = X*Y;
clustering = 1+floor((k*rand(N, 1)));    
end

