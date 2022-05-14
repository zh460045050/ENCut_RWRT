//
//  optimal_b.cpp
//  NCutH
//
//  Created by 朱磊 on 2018/8/5.
//  Copyright © 2018年 朱磊. All rights reserved.
//

#include "optimal_b.hpp"

optimalB::optimalB(){};

optimalB::optimalB(vector<double> data, double sigma)
{
    this->data = data;
    this->sigma = sigma;
    this->num_data = data.size();
    this->NCut = vector<double>(num_data, 0);
    this->ass_l = vector<double>(num_data, 0);
    this->ass_r = vector<double>(num_data, 0);
    this->cut = vector<double>(num_data, 0);
}

vector<double> optimalB::calEG(vector<double> X)
{
    vector<double> EG = vector<double>(num_data, 0);
    double X_mid = X[num_data / 2];
    for(int i=0; i<X.size(); i++)
        EG[i] = exp((X[0] - X[i]) / sigma);
    return EG;
}

vector<double> optimalB::calCS(vector<double> X)
{
    vector<double> CS = vector<double>(num_data, 0);
    double sum = 0;
    for(int i=0; i<X.size(); i++)
    {
        sum = sum + X[i];
        CS[i] = sum;
    }
    return CS;
}

vector<double> optimalB::calInv(vector<double> X)
{
    vector<double> Inv = vector<double>(num_data, 0);
    for(int i=0; i<X.size(); i++)
        Inv[i] = 1.0 / double(X[i]);
    return Inv;
}

void optimalB::calNCut()
{
    int n = data.size() - 1;
    vector<double> EG = calEG(data);
    vector<double> EGinv = calInv(EG);
    vector<double> CS = calCS(EG);
    vector<double> CSinv = calCS(EGinv);
    vector<double> EG_CSinv = vector<double>(num_data, 0);
    vector<double> EGinv_CS_in = vector<double>(num_data, 0);
    for(int i=0; i<EG.size(); i++)
    {
        EG_CSinv[i] = EG[i] * CSinv[i];
        EGinv_CS_in[i] = EGinv[i] * (CS[n] - CS[i]);
    }
    vector<double> CS_CSinv = calCS(EG_CSinv);
    vector<double> CSinv_CS_in = calCS(EGinv_CS_in);
    for(int k=0; k<data.size(); k++)
    {
        cut[k] = CSinv[k] * (CS[n] - CS[k]);
        ass_l[k] = CSinv_CS_in[k] + CSinv[k];
        ass_r[k] = CS_CSinv[n] - CS_CSinv[k] + CSinv_CS_in[n] - CSinv_CS_in[k];
        NCut[k] = cut[k] * (1.0 / ass_l[k] + 1.0 / ass_r[k]);
    }
}

void optimalB::calOptimal()
{
    std::vector<double>::iterator smallest = std::min_element(NCut.begin(), NCut.end());
    k = smallest - NCut.begin();
    b = (data[k] + data[k + 1]) / 2.0;
    min_ncut = double(*smallest);
}
