function [ clustering ] = initialKMeans( img, k )
%INITIALNCUT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[X, Y, Z] = size(img);
img_val = double(reshape(img, X*Y, Z));
clustering = kmeans(img_val,k);

end

