#pragma once

using scalar = double;
using label  = int;

#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <cassert>
#include <set>
#include <fstream>
// #include <boost/multi_array.hpp>
#include <unsupported/Eigen/CXX11/Tensor>

#include "Vector.hpp"


using svector = Vector<scalar>;
using Point = svector;

// using array_4 = boost::multi_array<scalar, 4>;
// using array_3 = boost::multi_array<scalar, 3>;
// using array_2 = boost::multi_array<scalar, 2>;

using array_4 = Eigen::Tensor<scalar, 4>;
using array_3 = Eigen::Tensor<scalar, 3>;
using array_2 = Eigen::Tensor<scalar, 2>;

#ifndef _SIZE_T
typedef __SIZE_TYPE__ label;
#endif


struct Parameter
{
    scalar x0, y0, x1, y1;  // Domain
    label nx, ny;           // Number of cells

    label k;                // Order of the scheme
    label dimPk;            // Dimension of the coefficient field
    label numGLP;           // Number of Gauss-Lobatto points

    scalar gamma;           // Heat capacity ratio

};


void readPara(std::string paraDir, Parameter& param);