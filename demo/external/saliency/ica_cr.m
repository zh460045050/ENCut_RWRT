function [yy]=ica_cr(im,r)
% color reference ica-r
% im          image
% r           reference signal
[h,w,d]=size(im);
im=im2double(im);
r=im2double(r);
len=h*w;
x=zeros(d,len);
for i=1:d,
    im1=im(:,:,i);
    x(i,:)=im1(:);
end;
[ptx,pty]=find(r>0);
[h,w,d]=size(r);
for i=1:d,
    im1=r(:,:,i);
    r1(i,:)=im1(:);
end;
r=r1;
 % remove the mean center the mixture  去均值
    x = x - mean( x, 2 )* ones( 1, len );
    % whiten the observed signals  白化观测信号 使得其自相关（协方差矩阵）为单位阵
    [ Ex, Dx ] = eig( cov( x' ) );%EVD decomposition [eigenvalues,eigenvectors]
    Q = Ex*sqrt( inv( Dx ) ) * Ex';%modified by li
    x = Q * x;
    [ m, n ] = size( x );
    %generate gaussian random variety v    产生高斯随机向量
    for i=1:d
    v = randn( 1, len );
    v = v - mean( v ) * ones( size( v ) );
    v = v / sqrt( var ( v ) );
    v1(i,:)=v;
    end;
    v=v1;

    % seting and initializing some parameters   设置初始参数
    gamma = 0.3;
    w_new = zeros(m,d);
    w_old = zeros(m, d);
    w_plus = zeros(m, d);
    mu = 0.4;
    max_iterat = 100;    %最大迭代次数
    iterat = 0;          %初始值
    rho = 0.7;
    xi = 0.001;          %门限
    eta = 0.5;
    y = zeros( d, len ); %输出信号
    Delta_L = zeros( m, 1);
    epsilon_IPI = 2e-2;   %相应判决门限
    look = [];
    ort = 1;
    % 主要迭代过程。
    while iterat < max_iterat && abs(ort-1)<1e-10
        y = w_new' * x;%signal object
        G_y = exp( -y.^2 / 2 );
        G_v = exp( -v.^2 / 2 );%reference gaussian
        rho =  mean( G_y,2) - mean( G_v ,2);
        g_w = (y-r)*(y-r)' - xi;
        mu = max( 0, mu + gamma * g_w );%
        %wk+1=wk-...
        D_G   =  -y .* exp( -y.^2 / 2 );%first derivative of gaussian
        %second derivative of gaussion
        D_D_G = ( y.^2 - ones( size( y ) ) ) .* exp( -y .^2 / 2 );
        D_g   =  2*(y-r);%first derivative of g(w)
        D_D_g =  2;%second derivative of g(w)
        %L'wk=...
        D_L   =  rho * mean( x .*( ones( m, 1) * D_G), 2) - 0.5 * mu  * mean( x.*(ones(m,1)*D_g), 2);
        D_w   =  rho * mean( D_D_G ) - 0.5 * mu * D_D_g;
        w_old = w_new;
        w_plus = w_old -eta*D_L/D_w;
        w_new = w_plus./norm(w_plus,2);
        iterat = iterat + 1;
        ort = w_old'*w_new;
    end
    y = w_new'*x;
    yy=(reshape(y,h,w));
    
%     level=graythresh(yy)*255+10;
%     yy=GuiYiHua(yy);
%     icabw=((yy>level)); 
%     
%     se = strel('square',5);
%     icabw=imerode(icabw,se);
%     icabw=bwselect(icabw,pty,ptx);
%     
%     icabw=imdilate(icabw,se);  
    
    
   