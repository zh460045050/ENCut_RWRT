//
//  mexFunction.cpp
//  gradDRW
//
//  Created by ?±ç? on 2018/3/25.
//  Copyright Â© 2018å¹??±ç?. All rights reserved.
//

#include <iostream>
#include "mex.h"
#include "bfgs.hpp"
using namespace std;
using namespace cv;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//读取输入数据
	double *inV = (double*)mxGetPr(prhs[0]);//矩阵维度
    double *inX = (double*)mxGetPr(prhs[1]);//矩阵维度
    double max_iter = mxGetScalar( prhs[2] );//最大迭代次数
    double lamda = mxGetScalar( prhs[3] );//步长
	double sigma = mxGetScalar( prhs[4] );//高斯参数

    int V_rows = mxGetM(prhs[0]); 
    int V_cols = mxGetN(prhs[0]);
	int X_rows = mxGetM(prhs[1]); 
	int X_cols = mxGetN(prhs[1]);

    int N = X_rows;
    int d = X_cols;

    Mat X(N, d, CV_64F);
    Mat V(d, 1, CV_64F);

    for(int i=0; i<V_rows; i++)
        for(int j=0; j<V_cols; j++)
            V.at<double>(i, j) = *(inV + i + j * V_rows);


    for(int i=0; i<X_rows; i++)
        for(int j=0; j<X_cols; j++)
            X.at<double>(i, j) = *(inX + i + j * X_rows);


    bfgs opt(V, X, max_iter, lamda, sigma);
    Mat optV = opt.optV;
    double optb = opt.optb;
    double min_ncut = opt.ncut;



    plhs[0] = mxCreateDoubleMatrix(optV.rows, optV.cols, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

    double *dV = mxGetPr(plhs[0]);
    double *b = mxGetPr(plhs[1]);
    double *ncut = mxGetPr(plhs[2]);
    *b = optb;
    *ncut = min_ncut;
    
    for (int i = 0; i < optV.rows; i++)
    {
        for (int j = 0; j < optV.cols; j++)
        {
            *(dV + i + j * optV.rows) = (double)optV.at<double>(i, j);
        }
    }
    
}
