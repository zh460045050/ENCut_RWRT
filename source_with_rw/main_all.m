%Example of normalized cut (NC) or average association (AA) clustering using our kernel cut 
% % The objective is NC or AA ONLY. No MRF term here so just use unary bounds.
% generate ring data, kernel and initial clustering
clear all;
close all;
clc;

energy_type = 'NCRW'; % 'NC' or 'RW' or 'NCRW'
weight_type = 'ENCut'; % 'ENCut' or 'NCut'
init_type = 'NCut'; % 'ENCut' or 'NCut' or 'Kmeans'
k = 25;
sigma = 50;
maxiter = 5;


img = imread('/Users/zhulei/DataSet/Segmentation/pictures/41004.jpg');
%img = imread('2007_003330.jpg');
%img = imread('3096.jpg');
img = imresize(img, 1.0);
[X, Y, Z] = size(img);
%w=fspecial('gaussian',[3 3], 1);
%img=imfilter(img,w);

%初始化相似矩阵
A = formWeight(img, sigma, 'ENCut');

%初始化标签
[clustering] = initialNCut(A, k);
%[clustering] = initialRandom(img, k);
%[clustering] = initialKMeans(img, k);

%显示初始标签
sc_labels = reshape(clustering, [X, Y]);
showFigs(img, sc_labels);

[~, labels] = moveMaking(A, k, clustering, maxiter, energy_type);

labels = reshape(labels, [X, Y]);
showFigs(img, labels);


