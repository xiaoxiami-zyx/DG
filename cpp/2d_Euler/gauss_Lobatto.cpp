#include "gauss_Lobatto.hpp"

Gauss_Lobatto::Gauss_Lobatto(label numGLP)
{
    numGLP_ = numGLP;
    lambda_.resize(numGLP_);
    weight_.resize(numGLP_);


    if(numGLP_ == 5)
    {
        lambda_[1] = -0.9061798459386639927976269;     
        lambda_[2] = -0.5384693101056830910363144;     
        lambda_[3] = 0;                                 
        lambda_[4] = 0.5384693101056830910363144;     
        lambda_[5] = 0.9061798459386639927976269;     
    
        weight_[1] = 0.2369268850561890875142640;
        weight_[2] = 0.4786286704993664680412915;
        weight_[3] = 0.5688888888888888888888889;
        weight_[4] = 0.4786286704993664680412915;
        weight_[5] = 0.2369268850561890875142640;
        
        //  Gauss-Lobatto Points
        lambdaL_[1] = -1;
        lambdaL_[2] = -0.6546536707079771437983;
        lambdaL_[3] = 0;
        lambdaL_[4] = 0.654653670707977143798;
        lambdaL_[5] = 1;
    }
}