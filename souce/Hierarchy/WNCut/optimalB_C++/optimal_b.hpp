//
//  optimal_b.hpp
//  NCutH
//
//  Created by 朱磊 on 2018/8/5.
//  Copyright © 2018年 朱磊. All rights reserved.
//

#ifndef optimal_b_hpp

#include <iostream>
#include <vector>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
using namespace std;
using namespace cv;

#define optimal_b_hpp


class optimalB
{
public:
    vector<double> data;
    double sigma;
    vector<double> NCut;
    vector<double> cut;
    vector<double> ass_l;
    vector<double> ass_r;
    double b;
    int k;
    int num_data;
    double min_ncut;
    optimalB();
    optimalB(vector<double> data, double sigma);
    vector<double> calEG(vector<double> X);
    vector<double> calCS(vector<double> X);
    vector<double> calInv(vector<double> X);
    void calNCut();
    void calOptimal();
};


#endif /* optimal_b_hpp */


