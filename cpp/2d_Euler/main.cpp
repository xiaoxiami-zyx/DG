#include "mesh.hpp"
#include "field.hpp"
#include "solver.hpp"

int main()
{
    Parameter param;
    readPara("var", param);

    Mesh mesh(param.x0, param.y0, param.x1, param.y1, param.nx, param.ny);

    // fieldCreator fieldCreator(mesh, param);

    // field u = fieldCreator.createField();
    
    // auto [nx, ny, dimPk] = u.dataShape();
    // auto [nx1, ny1, numGLP] = u.phiShape();
    // std::cout << "nx = " << nx << ", ny = " << ny << ", dimPk = " << dimPk << std::endl;
    // std::cout << "nx1 = " << nx1 << ", ny1 = " << ny1 << ", numGLP = " << numGLP << std::endl;

    Solver solver(mesh, param);
    solver.init();

    return 0;
}