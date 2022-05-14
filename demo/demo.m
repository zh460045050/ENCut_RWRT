clear;
clc;
addpath(genpath('encut_encode'));
addpath(genpath('experiment'));
addpath(genpath('external'));
%fileName = '253036.jpg';
%fileDir = 'img/';
img = imread('img/107793.jpg'); %image
k = 25; %number of partitions
thr = 0.3; %Threhod of the UCM

disp('This is the demo of the ENCut......')

%%%%%%Cluster-based ENCut Segmentation%%%%%%
clusterBased(img, k);

%%%%%%Cluster-based f-ENCut(fast exploring method) Segmentation%%%%%%
f_clusterBased(img, k);

%%%%%%Bipart-based ENCut Segmentation%%%%%%
bipartBased(img, k);

%%%%%%Hierichical-based ENCut Segmentation(with gPB)%%%%%%
hierBased(img, thr);

