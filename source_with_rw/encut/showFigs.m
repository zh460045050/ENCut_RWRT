function [ ] = showFigs( img, labels )
%SHOWRESULTS 此处显示有关此函数的摘要
%   此处显示详细说明

bd_result = show_result(img, labels);

imgSeg = colorSegmentedImage(img,labels);
figure; 
subplot(2,2,1); imagesc(img);
title('Source');
subplot(2,2,2); imagesc(bd_result);
title('Boundary');
subplot(2,2,3); imagesc(labels);
title('Labels');
subplot(2,2,4); imagesc(imgSeg);
title('Colorized');

end

