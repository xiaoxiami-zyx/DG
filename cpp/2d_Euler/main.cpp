#include "mesh.hpp"

int main()
{
    Parameter param;
    readPara("var", param);

    Mesh mesh(param.x0, param.y0, param.x1, param.y1, param.nx, param.ny);

    

    return 0;
}