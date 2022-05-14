//
//  optimal_v.cpp
//  NCutH
//
//  Created by 朱磊 on 2018/8/6.
//  Copyright © 2018年 朱磊. All rights reserved.
//

#include "optimal_v.hpp"

using namespace std;
using namespace cv;


optimalV::optimalV(){}

optimalV::optimalV(Mat V, Mat X, double sigma)
{
    this->V_M = V;
    this->X_M = X;
    this->sigma = sigma;
    Mat norm_V = V.clone();
    for(int i=0; i<norm_V.rows; i++)
        for(int j=0; j<norm_V.cols; j++)
            norm_V.at<double>(i, j) = V.at<double>(i, j) / norm(V, 2);
    Mat VX_M = X * norm_V;
    for(int i=0; i<VX_M.rows; i++)
        sortX.push_back(VX_M.at<double>(i, 0));
    sort(sortX.begin(), sortX.end());
    optb = optimalB(sortX, sigma);
}

void optimalV::calOptimalB()
{
    optb.calNCut();
    optb.calOptimal();
}

Mat optimalV::caldVX()
{
    Mat dVX;
    dVX = 1 / norm(V_M, 2) * X_M - 1 / pow(norm(V_M, 2), 3) * X_M * (V_M * V_M.t());
    return dVX;
}

Mat optimalV::caldNCut()
{
    int n = sortX.size() - 1;
    Mat dNCut_M = Mat::zeros(1, n+1, CV_64F);
    
    int k = optb.k;
    vector<double> E = calE(sortX, k);
    vector<double> Einv = calInv(E);
    vector<double> CS = calCS(E);
    vector<double> CSinv = calCS(Einv);
    
    vector<double> dNCut = vector<double>(sortX.size(), 0);
    vector<double> dCut = vector<double>(sortX.size(), 0);
    vector<double> dAssl = vector<double>(sortX.size(), 0);
    vector<double> dAssr = vector<double>(sortX.size(), 0);
    
    for(int l=0; l<=k; l++)
    {
        dCut[l] = 1.0 / sigma * E[l] * (CSinv[n] - CSinv[k]);
        
        if(l == 0)
            dAssl[l] = 0;
        else
            dAssl[l] = -2.0 * sigma * CS[l-1] / E[l];
        dAssl[l] += 2.0 / sigma * E[l] * (CSinv[k] - CSinv[l]);
        dAssl[l] += 1.0 / sigma * E[l] * (CSinv[n] - CSinv[k]);
        
        dAssr[l] = 1.0 / sigma * E[l] * (CSinv[n] - CSinv[k]);
        
        dNCut[l] = dCut[l] * (1.0 / optb.ass_l[k] + 1.0 / optb.ass_r[k]);
        dNCut[l] -= optb.cut[k] * (dAssl[l] / pow(optb.ass_l[k], 2) + dAssr[l] / pow(optb.ass_r[k], 2));
        dNCut_M.at<double>(0, l) = (double)dNCut[l];
    }
    
    for(int l=k+1; l<=n; l++)
    {
        dCut[l] = -1.0 / sigma * CS[k] / E[l];
        
        dAssl[l] = -1.0 / sigma * CS[k] / E[l];
        
        dAssr[l] = 2.0 / sigma * E[l] * (CSinv[n] - CSinv[l]);
        dAssr[l] -= 2.0 / sigma * (CS[l-1] - CS[k]) / E[l];
        dAssr[l] -= 1.0 / sigma * CS[k] / E[l];
        
        dNCut[l] = dCut[l] * (1.0 / optb.ass_l[k] + 1.0 / optb.ass_r[k]);
        dNCut[l] -= optb.cut[k] * (dAssl[l] / pow(optb.ass_l[k], 2) + dAssr[l] / pow(optb.ass_r[k], 2));
        dNCut_M.at<double>(0, l) = (double)dNCut[l];
    }

    return dNCut_M;
}

void optimalV::caldV()
{
    calOptimalB();
    Mat dVX = caldVX();
    Mat dNCut = caldNCut();
    dV = dNCut * dVX;
}

void bfgs()
{
    
}

vector<double> optimalV::calE(vector<double> X, int k)
{
    vector<double> E = vector<double>(X.size(), 0);
    for(int i=0; i<X.size(); i++)
        E[i] = 2 * (exp((X[i] - X[k])) / sigma) * (exp((X[i] - X[k])) / sigma);
    return E;
}

vector<double> optimalV::calCS(vector<double> X)
{
    vector<double> CS = vector<double>(X.size(), 0);
    double sum = 0;
    for(int i=0; i<X.size(); i++)
    {
        sum = sum + X[i];
        CS[i] = sum;
    }
    return CS;
}

vector<double> optimalV::calInv(vector<double> X)
{
    vector<double> Inv = vector<double>(X.size(), 0);
    for(int i=0; i<X.size(); i++)
        Inv[i] = 1.0 / double(X[i]);
    return Inv;
}

