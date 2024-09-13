#pragma once
#include "param.hpp"
#include "field.hpp"


// 将物理场转换为基函数的系数场
class Phi2Coff
{
    const std::vector<scalar>&   weight_;    
    const array_3&     phiGauss_;
    const std::vector<scalar>& mm_;


public:
    Phi2Coff(const std::vector<scalar>& weight, const array_3& phiGauss, const std::vector<scalar>& mm_);

    ~Phi2Coff(){};


    void computeCoff(Field& a);

};