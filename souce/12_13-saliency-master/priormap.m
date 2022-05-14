function sal = priormap(lab_im, ind, superlabel,pb_im,rgb_im)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Input parameter:
%     lab_im    : lab color image
%     ind       : positions of pixels inside the hull
%     superlabel: superpixel label
%     pb_im     : pb boundary image
%     rgb_im    : rgb color image
%Output parameter:
%     sal       : saliency map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lamda = 1; % the coefficient banlance the color and distance
theta = 1;
L = lab_im(:,:,1);
A = lab_im(:,:,2);
B =lab_im(:,:,3); % 3 color channels
[row,col] = size(L);
out_ind = setdiff(1:row*col,ind);
out_mean = [mean(L(out_ind)),mean(A(out_ind)),mean(B(out_ind))];
 
%% check whether the super pixel is completely inside convext

STATS = regionprops(superlabel, 'all');
sup_num = numel(STATS);
innersuper = [];
insup_mean = [];
for r = 1:sup_num %% check the superpixel not along the image sides
    indxy = STATS(r).PixelIdxList;
    if (numel(intersect(indxy,ind)) > 0.4 * numel(indxy))
        innersuper = [innersuper,r];
        insup_mean = [insup_mean;mean(L(indxy)),mean(A(indxy)), mean(B(indxy))];%%color of inner superpixel
    end   
end
outersuper = setdiff(1:sup_num,innersuper);%
outsup_mean = [];
for r = 1:numel(outersuper)
    out_indxy = STATS(outersuper(r)).PixelIdxList;
    outsup_mean = [outsup_mean;mean(L(out_indxy)),mean(A(out_indxy)), mean(B(out_indxy))];
end
innersup_num  = numel(innersuper);
pos_mat(sup_num,2) = 0;
color_mat(sup_num,3) = 0;
boundary_weight(innersup_num)=0;%%boundary weight map
pb_weights(innersup_num) = 0;
%%obtain boundary weight map
for i = 1:innersup_num
    mask= superlabel==innersuper(i);
    boundary = edge(mask).*pb_im;
    boundary_idx = find(boundary);
    pb_weights(i) = sum(sum(boundary))/sum(sum(edge(mask)));%%边缘pb值的均值
    boundary_weight(i) = sum((insup_mean(i,:) - out_mean).*(insup_mean(i,:) - out_mean))*pb_weights(i);
end
for m = 1:sup_num
    pixelind = STATS(m).PixelIdxList;
    indxy = STATS(m).PixelList;% indxy有两列 第一列是列数，第二列是行数
    pos_mat(m,:) = [mean(indxy(:,1)),mean(indxy(:,2))];
    color_mat(m,:) = [mean(L(pixelind)),mean(A(pixelind)),mean(B(pixelind))];
end   
%% compute color distance and spacial distance for prior map
pos_mat = (pos_mat - min(pos_mat(:)))/(max(pos_mat(:)) - min(pos_mat(:)));
color_mat = (color_mat - min(color_mat(:)))/(max(color_mat(:)) - min(color_mat(:)));

%%Compute color
mat_temp(sup_num,innersup_num) = 0;
vector_temp(sup_num) = 0;
for n = 1:innersup_num
    harris_sp_label = innersuper(n);
    harris_sp_color = color_mat(harris_sp_label,:);
    harris_sp_pos = pos_mat(harris_sp_label,:);    
    for q = 1:sup_num
        cur_sp_color = color_mat(q,:);
        cur_sp_pos = pos_mat(q,:);
        if(harris_sp_label == q)
            sal_temp = 0;
        else
            d_color = sqrt(sum((harris_sp_color - cur_sp_color).^2 ));
            d_space = sqrt(sum((harris_sp_pos - cur_sp_pos).^2));
            sal_temp = boundary_weight(n)/(d_color + lamda * d_space);
        end
        vector_temp(q) = sal_temp;
    end
    mat_temp(:,n) = vector_temp;
end
for n = 1:innersup_num
    harris_sp_label = innersuper(n);
    mat_temp(harris_sp_label,harris_sp_label) = mean(mat_temp(harris_sp_label,:));
end

%% 对应到图像中的superpixel中
sal_vector = mean(mat_temp,2);
sal(row,col) = 0;
for m = 1:sup_num
    pixelind = STATS(m).PixelIdxList;
    sal(pixelind) = sal_vector(m);   
end 
sal = (sal - min(sal(:)))/(max(sal(:)) - min(sal(:)));