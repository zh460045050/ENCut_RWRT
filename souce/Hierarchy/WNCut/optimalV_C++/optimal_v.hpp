//
//  optimal_v.hpp
//  NCutH
//
//  Created by 朱磊 on 2018/8/6.
//  Copyright © 2018年 朱磊. All rights reserved.
//

#ifndef optimal_v_hpp

#include <iostream>
#include <vector>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "optimal_b.hpp"
using namespace std;
using namespace cv;

#define optimal_v_hpp

class optimalV
{
public:
    Mat V_M;
    Mat X_M;
    double sigma;
    
    vector<double> sortX;
    optimalB optb;
    Mat dV;
    
    int num_data;
    
    optimalV();
    optimalV(Mat V, Mat X, double sigma);
    
    Mat caldVX();
    void calOptimalB();
    Mat caldNCut();
    void caldV();
    vector<double> calE(vector<double> X, int k);
    vector<double> calCS(vector<double> X);
    vector<double> calInv(vector<double> X);
    void bfgs();
};

#endif /* optimal_v_hpp */
