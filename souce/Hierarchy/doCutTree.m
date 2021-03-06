function [ labels, costList ] = doCutTree(img, k, type)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
%G = fspecial('gaussian', [5 5], 2);img = imfilter(img,G,'same');

img = im2double(img);
r = img(:,:,1);
g = img(:,:,2);
b = img(:,:,3);
[X, Y, Z] = size(img);
sigma = 50;

life_node = 1;
count_lab = 0;
labels = zeros(X, Y);
queue = {};
Vlist = {};
costList = [];
count_V = 0;

preV = ones(1, 3);

if length(type) == 5 && sum(type == 'ENCut') == 5
    salmap = doSaliency(img);
end

if length(type) == 5 && sum(type == 'SNCut') == 5
    salmap = doSaliency(img);
end


lab = 0;

labels = zeros(X, Y);
if length(type) == 5 && sum(type == 'RWCut') == 5
    [ d, son_idx, ncut ] = biRWCut( img, 5, labels, lab );
elseif length(type) == 5 && sum(type == 'ENCut') == 5
    [ d, son_idx, ncut ] = biENCut( img, 5, labels, lab, salmap);
elseif length(type) == 5 && sum(type == 'SNCut') == 5
    [ d, son_idx, ncut ] = biSNCut( img, labels, lab, salmap);
elseif length(type) == 3 && sum(type == 'Cut') == 3
    [ d, son_idx, ncut ] = biCut( img, labels, lab);
elseif length(type) == 4 && sum(type == 'NCut') == 4
    [ d, son_idx, ncut ] = biNCut( img, labels, lab);
elseif length(type) == 6 && sum(type == 'CCBCut') == 6
    [ d, son_idx, ncut ] = biCCBCut( img, labels, lab);
end

queue{1, 1} = -ncut; 
queue{1, 2} = son_idx;
queue{1, 3} = lab;
queue{1, 4} = d;
count_que = 1;

disp(strcat('doing...', type));
while life_node < k && length(queue) ~= 0
    if queue{count_que, 1} == 100
        break;
    end
    ncut = -queue{count_que, 1};
    
    son_idx = queue{count_que, 2};
    
    
    lab = queue{count_que, 3};
    d = queue{count_que, 4};
    if length(queue{count_que}) == 5
        preV = queue{count_que, 5};
    end
    queue(count_que, :) = [];
    count_que = count_que - 1;
    
    if son_idx == -1
        continue;
    end
    
    life_node = life_node + 1;
    
    count_V = count_V + 1;
    Vlist{count_V} = preV;
    costList(count_V) = ncut;
   
    
    max_lab = max(max(labels));
    d(d < 0 ) = lab;
    d(d ~= lab ) = max_lab + 1;
    labels(son_idx) = d;
    
    if life_node == k
        break;
    end
%     
%     if ncut < 1e-4
%         continue;
%     end
%     
    lab_l = lab;
    lab_r = max_lab + 1;

    if length(type) == 5 && sum(type == 'RWCut') == 5
        [ d_l, son_idx_l, ncut_l ] = biRWCut( img, 5, labels, lab_l );
    elseif length(type) == 5 && sum(type == 'ENCut') == 5
        [ d_l, son_idx_l, ncut_l ] = biENCut( img, 5, labels, lab_l, salmap );
    elseif length(type) == 5 && sum(type == 'SNCut') == 5
        [ d_l, son_idx_l, ncut_l ] = biSNCut( img, labels, lab_l, salmap);
    elseif length(type) == 3 && sum(type == 'Cut') == 3
        [ d_l, son_idx_l, ncut_l ] = biCut( img, labels, lab_l );
    elseif length(type) == 4 && sum(type == 'NCut') == 4
        [ d_l, son_idx_l, ncut_l ] = biNCut( img, labels, lab_l ); 
    elseif length(type) == 6 && sum(type == 'CCBCut') == 6
        [ d_l, son_idx_l, ncut_l ] = biCCBCut( img, labels, lab_l );
    end

    count_que = count_que + 1;
    queue{count_que, 1} = -ncut_l; 
    queue{count_que, 2} = son_idx_l; 
    queue{count_que, 3} = lab_l; 
    queue{count_que, 4} = d_l;
    
    if length(type) == 5 && sum(type == 'RWCut') == 5
        [ d_r, son_idx_r, ncut_r ] = biRWCut( img, 5, labels, lab_r );
    elseif length(type) == 5 && sum(type == 'ENCut') == 5
        [ d_r, son_idx_r, ncut_r ] = biENCut( img, 5, labels, lab_r, salmap);
    elseif length(type) == 5 && sum(type == 'SNCut') == 5
        [ d_r, son_idx_r, ncut_r ] = biSNCut( img, labels, lab_r, salmap);
    elseif length(type) == 3 && sum(type == 'Cut') == 3
        [ d_r, son_idx_r, ncut_r ] = biCut( img, labels, lab_r );
    elseif length(type) == 4 && sum(type == 'NCut') == 4
        [ d_r, son_idx_r, ncut_r ] = biNCut( img, labels, lab_r );
    elseif length(type) == 6 && sum(type == 'CCBCut') == 6
        [ d_r, son_idx_r, ncut_r ] = biCCBCut( img, labels, lab_r );
    end
    count_que = count_que + 1;
    queue{count_que, 1} = -ncut_r; 
    queue{count_que, 2} = son_idx_r; 
    queue{count_que, 3} = lab_r; 
    queue{count_que, 4} = d_r;
    
    
    queue = sortrows(queue, 1);
end

labels = labels + 1;



end

