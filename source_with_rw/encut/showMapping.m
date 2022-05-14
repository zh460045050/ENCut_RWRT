function [ ] = showMapping( d, sz)
%SHOWMAPPING �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
h = figure();
set(h,'name', 'ENCut Embeddings', 'Numbertitle', 'off')
for i = 1:9
    rr = reshape(mapminmax(d(:,i), 0, 1), sz);
    subplot(3, 3, i); imshow(rr);
    title(strcat(num2str(i), 'th')); 
end

end

