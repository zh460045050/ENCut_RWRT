function B = ReadDAT(image_size,data_path)

%��ͼ��superpixel��DAT�ļ�,
%DAT��ʽ:ÿ�ĸ��ֽڴ洢һ������
%image_size  ͼ���С���� [m n] 
%data_path   DAT�ļ� 

row = image_size(1);
colomn = image_size(2);
fid = fopen(data_path,'r');
A = fread(fid, row * colomn, 'uint32')';
A = A + 1;
B = reshape(A,[colomn, row]);
B = B';
fclose(fid);