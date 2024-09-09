#pragma once

#include "param.hpp"

class Mesh
{
    scalar       x0_, y0_, x1_, y1_;
    size_t        nx_, ny_;
    scalar       dx_, dy_;

    std::vector<scalar>  x_, y_;   // x and y coordinates of the nodes
    std::vector<scalar>  xc_, yc_; // x and y coordinates of the cell centers

public:
    Mesh (scalar x0, scalar y0, scalar x1, scalar y1, size_t nx, size_t ny);

    ~Mesh(){};
        
    scalar x0() const { return x0_; }

    scalar y0() const { return y0_; }

    scalar x1() const { return x1_; }

    scalar y1() const { return y1_; }

    size_t nx() const { return nx_; }

    size_t ny() const { return ny_; }

    scalar dx() const { return dx_; }

    scalar dy() const { return dy_; }

    const std::vector<scalar>& x() const { return x_; }

    const std::vector<scalar>& y() const { return y_; }

    const std::vector<scalar>& xc() const { return xc_; }

    const std::vector<scalar>& yc() const { return yc_; }

};