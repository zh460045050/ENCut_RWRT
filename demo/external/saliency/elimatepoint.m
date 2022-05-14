function corner_im2 = elimatepoint(corner_im2,thresh)
% elimate the point the four edge of the image

elipos = [];
[row,col] = size(corner_im2);
[x,y] = ind2sub([row,col],find(corner_im2 == 1));
elipos_x  = find( ( x>(row - thresh) )...
                |( x < (thresh + 1))); 
elipos = [elipos ...
          (y(elipos_x)-1)*row + x(elipos_x)];
elipos_y  = find( (y > (col - thresh))...
                | (y < (thresh + 1)));
elipos = [elipos; ...
          (y(elipos_y)-1)*row + x(elipos_y)];
corner_im2 (elipos)=0; 

