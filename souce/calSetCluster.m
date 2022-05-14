function [  ] = calSetSeg(type_alg, type_set, k)
%CALSETBR 此处显示有关此函数的摘要
%计算在BSDS500数据集中的BR、UE及T
%   此处显示详细说明
%   SetBR : 数据集平均BR值
%   SetUE : 数据集平均UE值
%   SetT : 数据集平均运行时间
%   type : 算法种类


if length(type_set) == 4 && sum(type_set == 'MSRA') == 4
    dirOutput=dir(fullfile('dataSets/IMGs/MSRAimg/','*.jpg'));
end
if length(type_set) == 5 && sum(type_set == 'DRIVE') == 5
    dirOutput=dir(fullfile('dataSets/IMGs/DRIVEimg/', '*.tif'));
end
if length(type_set) == 5 && sum(type_set == 'STARE') == 5
    dirOutput=dir(fullfile('dataSets/IMGs/STAREimg/', '*.ppm'));
end
if length(type_set) == 4 && sum(type_set == 'BSDS') == 4
    dirOutput=dir(fullfile('dataSets/IMGs/BSDSimg/', '*.jpg'));
end
if length(type_set) == 8 && sum(type_set == 'BSDStest') == 8
    dirOutput=dir(fullfile('../data/BSDS500/images', '*.jpg'));
end


LengthFiles = length(dirOutput);
fileNames={dirOutput.name};


for i = 1:LengthFiles
%     if( i == 5)
%         break;
%     end
    disp(strcat(type_alg, '...', type_set, '...', num2str(i), '/', num2str(LengthFiles)));
    %count
    curname = fileNames{i};
    %curname
    A=isstrprop(curname,'digit');
    B=curname(A);
    C=str2num(B);
    curimg = loadSource(C, type_set);
    
    %curgt = loadGroundTruth(C, type_set);
    
    %%%Code%%%
    [rows, cols, ~] = size(curimg);
    %t1 = clock;
    
%     
%     if length(type_alg) == 6 && sum(type_alg == 'CCBCut') == 6
%         [curlabel] = doCCBCutTree(curimg, k);
%     else
%         [curlabel, ~] = doCutTree(curimg, k, type_alg);
%     end

    curlabel = doCluster(curimg, k, type_alg);

    write_dir = fullfile('../smallSetCluster/', type_alg);
    write_path = fullfile(write_dir, strcat(B, '.mat'));
    if ~exist(write_dir)
        mkdir(write_dir);
    end
    
    segs{1} = curlabel + 1;
    save(write_path, 'segs');
    %t2 = clock;

end

end

