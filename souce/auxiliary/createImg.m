function [ img ] = createImg()
%CREATEIMG �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

img = zeros(255, 255, 3);
r = zeros(255, 255, 3);
b = zeros(255, 255, 3);
for i = 1:255
    for j = 1:255
        img(i, j, 1) = i;
        img(i, j, 3) = 255- j;
        r(i, j, 1) = i;
        b(i, j, 3) = 255-j;
    end
end
img = uint8(img);
r = uint8(r);
b = uint8(b);
figure();imshow(r);
figure();imshow(b);
figure();imshow(img);

end

