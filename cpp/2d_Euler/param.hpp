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

#include "Vector.hpp"


using svector = Vector<scalar>;
using Point = svector;


struct Parameter
{
    scalar x0, y0, x1, y1;  // Domain
    label nx, ny;           // Number of cells

    label k;                // Order of the scheme
    label dimPk;            // Dimension of the coefficient field
    label numGLP;           // Number of Gauss-Lobatto points

};


void readPara(std::string paraDir, Parameter& param);