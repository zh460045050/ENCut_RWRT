
clear all
clc
addpath(genpath('./ColorFeatures'));

%imName = '0_18_18957.jpg';
imName = '/Users/zhulei/DataSet/Segmentation/pictures/3096.jpg';
lab_im = RGB2Lab(imread(imName));%%LAB image
rgb_im = imread(imName);%%RGB image
% pb_im = im2double(imread([imName(1:end-4) '_pb.png']));%% pb edge image after bi-segmenting by otsu
% sulabel_im = ReadDAT(size(pb_im),[imName(1:end - 4) '.dat']);%Superpixel Label


%% wseg from DRW
%%%%%%%%%%pb%%%%%%%%%%%%%%%
pb = pbCGTG_nonmax(double(rgb_im)/255); %superpixel segmentation-->probability img
[f_im1, ~] = max(pb, [], 3);
figure,imshow(f_im1);
pb = f_im1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = watershed(pb); 
wr = L == 0;
aa = 1 - wr;
[labeled,numObjects] = bwlabel(aa, 8);
STATS = regionprops(labeled,'centroid');
centroids = cat(1, STATS.Centroid);

k = max(max(labeled));
p1 = round(centroids(:,1));
p0 = round(centroids(:,2));

%%%%DRW+water+pb%%%%%
%[sulabel_im, b] = mex_DRW(rgb_im, 200, 1e-5, f_im1, 5, p0, p1, k);
%[sulabel_im, b, pb_seeds_x, pb_seeds_y] = pb_seed(rgb_im, 2000, 1e-4, f_im1, 0);
[sulabel_im, area_seeds_x, area_seeds_y] = mexDRW(rgb_im, 1500, 1e-3, 50);

%[sulabel_im, numlabels] = slicmex(rgb_im, 500, 20);

%[sulabel_im] = mex_ers(double(rgb_im),500,0.5,5.0,1);

length(unique(sulabel_im))

result = show_result(rgb_im, sulabel_im);
figure,imshow(result);


%% Obtain interesting points
thresh = 26; % for elimate the side point
corner_im2 = getsalientpoints(rgb_im);
corner_im = elimatepoint(corner_im2,thresh); % elimate the points closing to the boundary of images
%% Calculate prior map
[row,col] = size(corner_im);
[y,x] = ind2sub([row,col],find(corner_im == 1));
dt = DelaunayTri(x,y);
if(~size(dt,1))
    return;
end
[k,av] = convexHull(dt);
BW = roipoly(corner_im,x(k),y(k));
pixel = regionprops(BW,'all');
ind = pixel.PixelIdxList;
out_ind = setdiff(1:row*col,ind);
sal_super = priormap(lab_im, ind, sulabel_im,pb,rgb_im);

%% Revised observation likelihood
convexhull = ReviseConvexHull(sal_super, rgb_im, sulabel_im);
hull_pixel = regionprops(convexhull,'all');
ind = [];
out_ind = [];
for i = 1:length(hull_pixel)
   ind = [ind;hull_pixel(i).PixelIdxList];%%�ҵ�mask�İ�ɫ����
end
out_ind = setdiff(1:row*col,ind);
[PrI_sal, PrI_bk,PrO_sal,PrO_bk] = likelihoodprob(rgb_im, ind,out_ind);

%% Bayesian combination
psal_I = sal_super(ind);
psal_O = sal_super(out_ind');
Pr_0=(PrI_sal.*psal_I)./(PrI_sal.*psal_I+PrI_bk.*(1 - psal_I));%so called saliency ���ڵ�saliency
Pr_B=(PrO_sal.*psal_O)./(PrO_sal.*psal_O+PrO_bk.*(1-psal_O));%so called saliency �����?aliency
saliencymap = zeros(row,col);
saliencymap(ind) = Pr_0;
saliencymap(out_ind) = Pr_B;
saliencymap = (saliencymap - min(saliencymap(:)))/(max(saliencymap(:)) - min(saliencymap(:)));
figure,imshow(saliencymap);
