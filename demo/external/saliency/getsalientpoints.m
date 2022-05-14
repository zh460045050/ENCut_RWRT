function corner_im2 = getsalientpoints(rgb_im)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%compute the harris point
%Input:  
%     rgb_im    : RGB color image
%Output: 
%     corner_im2: detected harris points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_g=1.5; % parameters for computing harris points
sigma_a=5;  % parameters for computing harris points
nPoints=30; % number of salient points 
rgb_im=double(rgb_im); 
Mboost = BoostMatrix(rgb_im);
boost_im= BoostImage(rgb_im,Mboost);
[EnIm]= ColorHarris(boost_im,sigma_g,sigma_a,0.04,1);
[x_max,y_max,corner_im2,num_max]=getmaxpoints(EnIm,nPoints);