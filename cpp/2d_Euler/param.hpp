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
#include <boost/multi_array.hpp>

#include "Vector.hpp"


using svector = Vector<scalar>;
using Point = svector;

using array_4 = boost::multi_array<scalar, 4>;
using array_3 = boost::multi_array<scalar, 3>;
using array_2 = boost::multi_array<scalar, 2>;

#ifndef _SIZE_T
typedef __SIZE_TYPE__ size_t;
#endif


struct Parameter
{
    scalar x0, y0, x1, y1;  // Domain
    size_t nx, ny;           // Number of cells

    size_t k;                // Order of the scheme
    size_t dimPk;            // Dimension of the coefficient field
    size_t numGLP;           // Number of Gauss-Lobatto points

};


void readPara(std::string paraDir, Parameter& param);