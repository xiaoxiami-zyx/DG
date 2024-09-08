#pragma once
#include "param.hpp"

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