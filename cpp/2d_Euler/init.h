#ifndef __INIT_H__
#define __INIT_H__

#include "param.hpp"


int region(scalar x, scalar y)
{
    if (x >= 0.8 && x <= 1.0 && y >= 0.8 && y <= 1.0)
        return 1;
    else if (x >= 0.0 && x <= 0.8 && y >= 0.8 && y <= 1.0)
        return 2;
    else if (x >= 0.0 && x <= 0.8 && y >= 0.0 && y <= 0.8)
        return 3;
    else if (x >= 0.8 && x <= 1.0 && y >= 0.0 && y <= 0.8)
        return 4;

    assert(false);
}


scalar rho0(scalar x, scalar y)
{
    int r = region(x, y);
    if (r == 1)
        return 1.5;
    else if (r == 2)
        return 0.5323;
    else if (r == 3)
        return 0.138;
    else if (r == 4)
        return 0.5323;
}

scalar u0(scalar x, scalar y)
{
    int r = region(x, y);
    if (r == 1)
        return 0.0;
    else if (r == 2)
        return 1.206;
    else if (r == 3)
        return 1.206;
    else if (r == 4)
        return 0.0;
}

scalar v0(scalar x, scalar y)
{
    int r = region(x, y);
    if (r == 1)
        return 0.0;
    else if (r == 2)
        return 0.0;
    else if (r == 3)
        return 1.206;
    else if (r == 4)
        return 1.206;
}

scalar p0(scalar x, scalar y)
{
    int r = region(x, y);
    if (r == 1)
        return 1.5;
    else if (r == 2)
        return 0.3;
    else if (r == 3)
        return 0.029;
    else if (r == 4)
        return 0.3;
}


#endif