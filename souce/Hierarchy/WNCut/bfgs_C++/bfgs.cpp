//
//  bfgs.cpp
//  NCutH
//
//  Created by 朱磊 on 2018/8/6.
//  Copyright © 2018年 朱磊. All rights reserved.
//

#include "bfgs.hpp"
bfgs::bfgs(){}
bfgs::bfgs(Mat V, Mat X, int max_iter, double lamda, double sigma)
{
    this->V = V;
    this->X = X;
    this->max_iter = max_iter;
    this->lamda = lamda;
    this->sigma = sigma;
    dobfgs();
}
void bfgs::dobfgs()
{
    int N = X.rows;
    int d = X.cols;
    Mat D_i = Mat::eye(d, d, CV_64F);
    Mat I = Mat::eye(d, d, CV_64F);
    Mat V_i = V;
    optimalV optv(V_i, X, sigma);
    optv.caldV();
    Mat g_i = optv.dV.t();
    int iter = 0;
    Mat d_i;
    Mat y_i;
    Mat g_ni;
    Mat s_i;
    while(iter < max_iter)
    {
        d_i = -D_i * g_i;
        s_i = lamda * d_i;
        V_i = V_i + s_i;
        optimalV optv(V_i, X, sigma);
        optv.caldV();
        g_ni = optv.dV.t();
        y_i = g_ni - g_i;
        if(norm(g_i, 2) < 1e-5)
            break;
        D_i = (I -  (s_i * y_i.t()) / (y_i.t() * s_i)) * D_i * (I -  (y_i * s_i.t()) / (y_i.t() * s_i)) + (s_i * s_i.t()) / (y_i.t() * s_i);
        g_i = g_ni;
        iter ++;
    }
    optimalV final_optv(V_i, X, sigma);
    optv.caldV();
    this->optV = V_i.t();
    this->optb = optv.optb.b;
    this->k = optv.optb.k;
    this->ncut = optv.optb.min_ncut;
}
