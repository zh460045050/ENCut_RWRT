//
//  mexFunction.cpp
//  gradDRW
//
//  Created by ?±ç? on 2018/3/25.
//  Copyright Â© 2018å¹??±ç?. All rights reserved.
//

#include <iostream>
#include "mex.h"
#include "optimal_b.hpp"
using namespace std;
using namespace cv;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//读取输入数据
	double *inX = (double*)mxGetPr(prhs[0]);//矩阵维度
	double sigma = mxGetScalar( prhs[ 1 ] );//高斯参数

	int rows = mxGetM(prhs[0]); 
	int cols = mxGetN(prhs[0]);
	int N = rows * cols;

	vector<double> X;
	for(int i=0; i<N; i++)
		X.push_back(*(inX + i));

	optimalB opB(X, sigma);
    opB.calNCut();

    plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, N, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, N, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1, N, mxREAL);

    double *NCut = mxGetPr(plhs[0]);
    double *cut = mxGetPr(plhs[1]);
    double *ass_l = mxGetPr(plhs[2]);
    double *ass_r = mxGetPr(plhs[3]);           
    
    for(int i=0; i<N; i++)
    {
    	*(NCut + i) = opB.NCut[i];
    	*(cut + i) = opB.cut[i];
    	*(ass_l + i) = opB.ass_l[i];
    	*(ass_r + i) = opB.ass_r[i];
    }
}
