function [ ] = viewLabel( labels )
%VIEWLABEL �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��


for i = 1:max(max(labels))
    rr = zeros(size(labels));
    rr(labels == i) = 1;
    figure();imshow(rr);
end

end

