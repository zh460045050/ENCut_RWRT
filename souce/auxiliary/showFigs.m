function [ ] = showFigs( img, labels )
%SHOWRESULTS �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

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

