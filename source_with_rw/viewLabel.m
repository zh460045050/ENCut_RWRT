function [ ] = viewLabel( labels )
%VIEWLABEL 此处显示有关此函数的摘要
%   此处显示详细说明


for i = 1:max(max(labels))
    rr = zeros(size(labels));
    rr(labels == i) = 1;
    figure();imshow(rr);
end

end

