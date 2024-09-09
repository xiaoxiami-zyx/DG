#include "mesh.hpp"

Mesh::Mesh(scalar x0, scalar y0, scalar x1, scalar y1, size_t nx, size_t ny):
    x0_(x0), y0_(y0), x1_(x1), y1_(y1), nx_(nx), ny_(ny), dx_((x1-x0)/nx), dy_((y1-y0)/ny)
{
    assert(x0 < x1);
    assert(y0 < y1);
    assert(nx > 0);
    assert(ny > 0);

    x_.resize(nx+1);
    y_.resize(ny+1);
    xc_.resize(nx);
    yc_.resize(ny); 

    for (size_t i=0; i<=nx; ++i)
        x_[i] = x0 + i*dx_;
    
    for (size_t j=0; j<=ny; ++j)
        y_[j] = y0 + j*dy_;
    
    for (size_t i=0; i<nx; ++i)
        xc_[i] = 0.5*(x_[i] + x_[i+1]);
    
    for (size_t j=0; j<ny; ++j)
        yc_[j] = 0.5*(y_[j] + y_[j+1]);

}
