function [ clustering ] = initialRandom( img, k )
%INITIALNCUT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[X, Y, Z] = size(img);
N = X*Y;
clustering = 1+floor((k*rand(N, 1)));    
end

