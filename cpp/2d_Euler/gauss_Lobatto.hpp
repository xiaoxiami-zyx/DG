#pragma once
#include "param.hpp"

/**
 * @brief 2D Gauss-Lobatto 积分点和权重
 * 
 * 基函数 1, X, Y, X^2-1/3, XY, Y^2-1/3
 * 区间 [-1, 1]
 */
class Gauss_Lobatto
{
    label numGLP_;        // Number of Gauss-Lobatto points
    std::vector<scalar> lambda_;   // 
    std::vector<scalar> weight_;   // Gauss-Lobatto weights
    std::vector<scalar> lambdaL_;  // Gauss-Lobatto points

public:
    Gauss_Lobatto(label numGLP);

    ~Gauss_Lobatto(){};

    label numGLP() const { return numGLP_; }

    const std::vector<scalar>& lambda() const { return lambda_; }

    const std::vector<scalar>& weight() const { return weight_; }

    const std::vector<scalar>& lambdaL() const { return lambdaL_; }

};