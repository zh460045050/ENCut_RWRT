function [  ] = calSetSeg(type_alg, type_set, k)
%CALSETBR �˴���ʾ�йش˺�����ժҪ
%������BSDS500���ݼ��е�BR��UE��T
%   �˴���ʾ��ϸ˵��
%   SetBR : ���ݼ�ƽ��BRֵ
%   SetUE : ���ݼ�ƽ��UEֵ
%   SetT : ���ݼ�ƽ������ʱ��
%   type : �㷨����


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

