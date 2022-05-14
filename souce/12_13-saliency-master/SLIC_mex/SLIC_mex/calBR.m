function [ BR ] = calBR( resultimg, groundtruth, sigma)
%CALBR 此处显示有关此函数的摘要
%计算评价标准boundary adherence(BE)
%   此处显示详细说明
%   resultimg : 超像素分割结果
%   groundtruth : 人工标记图片分割结果
%   sigma: 像素搜索范围

[rows, cols, ~] = size(resultimg);
count = 0;
error = 0;
go_size = 0;
go = zeros(sigma * sigma, 2);
for i = -sigma : sigma
    for j = -sigma : sigma
        go_size = go_size + 1;
        go(go_size, 1) = i;
        go(go_size, 2) = j;
    end
end
for i = 1 : rows
    for j = 1 : cols
        if groundtruth(i, j) == 0
           continue;
        else
            count = count + 1;
            flag = 0;
            for k = 1 : go_size
                xp = i + go(k, 1);
                yp = j + go(k, 2);
                if xp > 0 && yp>0 && xp <=rows && yp <= cols && resultimg(xp, yp) == 255
                    flag = 1;
                    break;
                end 
            end
            if flag == 0
                error = error + 1;
            end
        end
    end
end

BR = 1 - error / count;

end

