function convexhull = ReviseConvexHull(sal_super, rgb_im, sulabel_im)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The function is to calculate the revised convex hull by soft-segmentation
%% Input parameter:
%%     sal_super : prior map
%%     rgb_im    : rgb color image
%%     sulabel_im: superpixel label
%% Output parameter:
%%     convexhull: revised convex hull by combination of ICA_R and KDE
%%ICA MAP
[row col] = size(sal_super);
thresh = graythresh(sal_super);
bi_im = sal_super>=thresh;
SE1 = strel('square',11);
bi_im2 = imerode(bi_im,SE1,'same');
salhull = ica_cr(rgb_im,bi_im2);
ica_map_original = mat2gray(salhull);
thresh = graythresh(ica_map_original);
ica_map = ica_map_original>=thresh; 
ica_map = double(ica_map);
ica_map = imfill(ica_map);

%%KDE MAP
rgb_im = double(rgb_im);
r = rgb_im(:,:,1);
g = rgb_im(:,:,2);
b = rgb_im(:,:,3);
r = r./(r+g+b+eps);
g = g./(r+g+b+eps);
s = (r+g+b)/3;%%ͼ���rgs�ռ�
STATS = regionprops(sulabel_im,'all');
sup_num = numel(STATS);
insup_mean = [];
sup_mean = [];
ind = [];
bi_im = bwlabel(bi_im);
ss = regionprops(bi_im,'all');
for i = 1:length(ss)
    ind = [ind;ss(i).PixelIdxList];
end
for i = 1:sup_num
    indxy = STATS(i).PixelIdxList;
    sup_mean = [sup_mean;mean(r(indxy)),mean(g(indxy)), mean(s(indxy))];
    if (numel(intersect(indxy,ind)) > 0.4 * numel(indxy))
        insup_mean = [insup_mean;mean(r(indxy)),mean(g(indxy)), mean(s(indxy))];
        % region mean of each region color in RGB space
    end   
end
stdr=std(insup_mean(:,1));
stdg=std(insup_mean(:,2));
stdb=std(insup_mean(:,3));
% F_t = ksdensity(sup_mean,'NumPoints',351);
[~,F_t(:,1)] = ksdensity(sup_mean(:,1),sup_mean(:,1));
[~,F_t(:,2)]= ksdensity(sup_mean(:,2),sup_mean(:,2));
[~,F_t(:,3)] = ksdensity(sup_mean(:,3),sup_mean(:,3));
% F_t(:,2) = ksdensity(sup_mean(:,2),'NumPoints',190);
% F_t(:,3) = ksdensity(sup_mean(:,3),'NumPoints',190);
% F_t=kde_color1(sup_mean',insup_mean',[stdr stdg stdb]);%data location and bandwidth
F_t = max(F_t,2);
full_F_t = zeros(row,col);
for i = 1:sup_num
    indxy = STATS(i).PixelIdxList;
    full_F_t(indxy) = F_t(i);
end
thresh = graythresh(full_F_t);
kde_map = full_F_t>=thresh;
%% Combine two maps
deta = 0.5;
original_convexhull = ica_map.*(1-exp(-kde_map/deta));
if ~isempty(find(original_convexhull > 0))
    thresh = graythresh(original_convexhull);
    convexhull = original_convexhull>=thresh;
end
convexhull = double(convexhull);
convexhull = imfill(convexhull);
convexhull = bwlabel(convexhull);