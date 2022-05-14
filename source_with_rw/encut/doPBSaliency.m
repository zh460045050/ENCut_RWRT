function [ saliencymap ] = doPBSaliency( rgb_im, pb )
%DOSALIENCY ´Ë´¦ÏÔÊ¾ÓÐ¹Ø´Ëº¯ÊýµÄÕªÒª
%   ´Ë´¦ÏÔÊ¾ÏêÏ¸ËµÃ÷
lab_im = RGB2Lab(rgb_im);
%pb = pbCGTG_nonmax(double(rgb_im)/255); %superpixel segmentation-->probability img
[f_im1, ~] = max(pb, [], 3);
%figure,imshow(f_im1);
pb = f_im1;

L = watershed(pb); 
wr = L == 0;
aa = 1 - wr;
[labeled,numObjects] = bwlabel(aa, 8);
STATS = regionprops(labeled,'centroid');
centroids = cat(1, STATS.Centroid);

k = max(max(labeled));
p1 = round(centroids(:,1));
p0 = round(centroids(:,2));
[sulabel_im, area_seeds_x, area_seeds_y] = mexDRW(rgb_im, 1500, 1e-3, 50);

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
   ind = [ind;hull_pixel(i).PixelIdxList];%%ï¿½Òµï¿½maskï¿½Ä°ï¿½É«ï¿½ï¿½ï¿½ï¿½
end
out_ind = setdiff(1:row*col,ind);
[PrI_sal, PrI_bk,PrO_sal,PrO_bk] = likelihoodprob(rgb_im, ind,out_ind);

%% Bayesian combination
psal_I = sal_super(ind);
psal_O = sal_super(out_ind');
Pr_0=(PrI_sal.*psal_I)./(PrI_sal.*psal_I+PrI_bk.*(1 - psal_I));%so called saliency ï¿½ï¿½ï¿½Úµï¿½saliency
Pr_B=(PrO_sal.*psal_O)./(PrO_sal.*psal_O+PrO_bk.*(1-psal_O));%so called saliency ï¿½ï¿½ï¿½ï¿½ï¿?aliency
saliencymap = zeros(row,col);
saliencymap(ind) = Pr_0;
saliencymap(out_ind) = Pr_B;
saliencymap = (saliencymap - min(saliencymap(:)))/(max(saliencymap(:)) - min(saliencymap(:)));

end

