function [ groundtruth ] = loadGroundTruth(num)
%LOADGROUNDTRUTH 此处显示有关此函数的摘要
%读取BSD数据集中的人工标记
%   此处显示详细说明
%   num:图片编号
%   groundtruth:人工标记结果
data = load('gd/' + string(num) + '.mat');
cellgt = data(1).groundTruth;
groundtruth = cellgt{1}.Boundaries;
end

