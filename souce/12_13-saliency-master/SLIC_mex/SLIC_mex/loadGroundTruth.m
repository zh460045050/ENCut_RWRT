function [ groundtruth ] = loadGroundTruth(num)
%LOADGROUNDTRUTH �˴���ʾ�йش˺�����ժҪ
%��ȡBSD���ݼ��е��˹����
%   �˴���ʾ��ϸ˵��
%   num:ͼƬ���
%   groundtruth:�˹���ǽ��
data = load('gd/' + string(num) + '.mat');
cellgt = data(1).groundTruth;
groundtruth = cellgt{1}.Boundaries;
end

