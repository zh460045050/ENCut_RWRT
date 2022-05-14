function sImg = colorSegmentedImage(img, seg)
% This function outputs a new image with each segmented region containing
% the average color of that region in the original image.
% inputs:   an image
%           a segmentation mask seg
% outputs:  new image with each segmented region containing the average
%           color of that region in the original image

% author: Nathan D. Cahill
% email: nathan.cahill@rit.edu
% date: 21 February 2018

% zero out seg on boundaries
seg([1 end],:) = 0;
seg(:,[1 end]) = 0;

% initialize sImg
sImg = zeros(size(img,1),size(img,2),3);

% calculate the average intensity value for each segment based on the mask
% seg and set the value of all pixels in that segment equal to that average
% value
sel = strel('disk',1);
for i = 1:max(seg(:))
    mask = imerode(seg==i,sel);
    if ~any(mask(:))
        continue;
    end
    tmp1 = img(:,:,1);
    tmp2 = img(:,:,2);
    tmp3 = img(:,:,3);
    tmp1 = tmp1(:);
    tmp2 = tmp2(:);
    tmp3 = tmp3(:);
    muR = mean(tmp1(mask));
    muG = mean(tmp2(mask));
    muB = mean(tmp3(mask));
    %mask = reshape(mask,size(img,1),size(img,2));
    sImg(:,:,1) = sImg(:,:,1)+mask.*muR;
    sImg(:,:,2) = sImg(:,:,2)+mask.*muG;
    sImg(:,:,3) = sImg(:,:,3)+mask.*muB;
end

% set class of sImg to be class of img
sImg = feval(class(img),sImg);

end

