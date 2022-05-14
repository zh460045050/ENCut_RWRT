
%img = imread('/Users/zhulei/DataSet/Segmentation/pictures/41004.jpg');
img = imread('test5.png');
%img = imresize(img, 0.5);
% w=fspecial('gaussian',[3 3], 1);
% img=imfilter(img,w);
tic;
labels = ENRW(img, 3, 20, 50);
toc;    