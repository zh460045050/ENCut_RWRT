clear;
tic;

dirOutput = dir(fullfile('img','*.jpg'));   % 提取路径 
fileNames = {dirOutput.name}';   % 获得符合条件文件名

img = imread('/Users/zhulei/Documents/study/C Work/nodeRW/nodeRW/pictures/253036.jpg');
%for i = 1: length(fileNames)
%fileNames = '253036.jpg';
%img = imread(fullfile('img', fileNames));
[label, numlabels] = slicmex(img, 300   );
%label = delHole(label);
supxl_num = length(unique(label))
result = show_result(img, label);
%idx = find('.'==fileNames);
%imname = fileNames(1:idx-1);
%gd = loadGroundTruth(imname);
%calBR(result, gd, 2)
figure(1); imshow(result);
%end