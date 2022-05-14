//
//  bfgs.hpp
//  NCutH
//
//  Created by 朱磊 on 2018/8/6.
//  Copyright © 2018年 朱磊. All rights reserved.
//

#ifndef bfgs_hpp
#include <iostream>
#include <vector>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "optimal_v.hpp"
#define bfgs_hpp

using namespace std;
using namespace cv;

#endif /* bfgs_hpp */

class bfgs
{
    Mat V;
    Mat X;
    int max_iter;
    double lamda;
    double sigma;
    
public:
    Mat optV;
    double optb;
    int k;
    double ncut;
    
    bfgs();
    bfgs(Mat V, Mat X, int max_iter, double lamda, double sigma);
    void dobfgs();
};
