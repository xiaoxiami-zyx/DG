#pragma once
#include "param.hpp"
#include "gauss_Lobatto.hpp"

/**
 * @brief 基函数
 * 
 * 基函数 1, X, Y, X^2-1/3, XY, Y^2-1/3
 * 区间 [-1, 1]
 */
class BaseFunction
{
    label     dimPk_;  // 多项式阶数
    label     numGLP_; // Gauss-Lobatto 积分点数

    std::vector<scalar> phiRU_;   // 基函数在点(1,1)的值
    std::vector<scalar> phiLU_;   // 基函数在点(-1,1)的值
    std::vector<scalar> phiRD_;   // 基函数在点(1,-1)的值
    std::vector<scalar> phiLD_;   // 基函数在点(-1,-1)的值
    
    std::vector<scalar> mm_;  // mass matrix

    array_3     phiGauss_;     // [numGLP, numGLP, dimPk]  基函数在Gauss积分点的值
    array_4     phiGaussLL_;   // [numGLP, numGLP, dimPk, 3]
    array_2     phiGaussR_;    // [numGLP, dimPk]  right
    array_2     phiGaussU_;    // [numGLP, dimPk]  up
    array_2     phiGaussD_;    // [numGLP, dimPk]  down
    array_2     phiGaussL_;    // [numGLP, dimPk]  left
    array_3     phiGauss_dx_;  // [numGLP, numGLP, dimPk] 基函数在Gauss积分点的导数
    array_3     phiGauss_dy_;  // [numGLP, numGLP, dimPk] 基函数在Gauss积分点的导数

public:
    BaseFunction(label dimPk, label numGLP, Gauss_Lobatto& glp, scalar dx, scalar dy);

    ~BaseFunction(){};

    const array_3& phiGauss() const {return phiGauss_;};

    const std::vector<scalar>& mm () const {return mm_;};


    // 定义基函数
    static scalar basefunc_0(scalar x, scalar y) { (void)x; (void)y; return 1.0; };
    static scalar basefunc_1(scalar x, scalar y) { (void)y; return x; };
    static scalar basefunc_2(scalar x, scalar y) { (void)x; return y; };
    static scalar basefunc_3(scalar x, scalar y) { (void)y; return x * x - 1.0 / 3.0; };
    static scalar basefunc_4(scalar x, scalar y) { return x * y; };
    static scalar basefunc_5(scalar x, scalar y) { (void)x; return y * y - 1.0 / 3.0; };

    static scalar basefunc_dx_0(scalar x, scalar y, scalar dx, scalar dy) { (void)x; (void)y; (void)dx; (void)dy; return 0.0; };
    static scalar basefunc_dx_1(scalar x, scalar y, scalar dx, scalar dy) { (void)x; (void)y; (void)dy; return 1.0 / dx; };
    static scalar basefunc_dx_2(scalar x, scalar y, scalar dx, scalar dy) { (void)x; (void)y; (void)dx; (void)dy; return 0.0; };
    static scalar basefunc_dx_3(scalar x, scalar y, scalar dx, scalar dy) { (void)y; (void)dy; return 2.0 * x / dx; };
    static scalar basefunc_dx_4(scalar x, scalar y, scalar dx, scalar dy) { (void)x; (void)dy; return y / dx; };
    static scalar basefunc_dx_5(scalar x, scalar y, scalar dx, scalar dy) { (void)x; (void)y; (void)dx; (void)dy; return 0.0; };

    static scalar basefunc_dy_0(scalar x, scalar y, scalar dx, scalar dy) { (void)x; (void)y; (void)dx; (void)dy; return 0.0; };
    static scalar basefunc_dy_1(scalar x, scalar y, scalar dx, scalar dy) { (void)x; (void)y; (void)dx; (void)dy; return 0.0; };
    static scalar basefunc_dy_2(scalar x, scalar y, scalar dx, scalar dy) { (void)x; (void)y; (void)dx; return 1.0 / dy; };
    static scalar basefunc_dy_3(scalar x, scalar y, scalar dx, scalar dy) { (void)x; (void)y; (void)dx; (void)dy; return 0.0; };
    static scalar basefunc_dy_4(scalar x, scalar y, scalar dx, scalar dy) { (void)y; (void)dx; return x / dy; };
    static scalar basefunc_dy_5(scalar x, scalar y, scalar dx, scalar dy) { (void)x; (void)dx; return 2.0 * y / dy; };

};