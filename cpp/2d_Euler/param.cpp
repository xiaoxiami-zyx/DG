#include "param.hpp"

void readPara(std::string paraDir, Parameter& param)
{
    std::string invStr;
    std::ifstream inputPara{paraDir.c_str(), std::ios::in};
    if(!inputPara)
    {
        std::cerr<<paraDir<<" could not be opened"<<std::endl;
        exit(EXIT_FAILURE);
    }


    inputPara >> invStr >> param.x0;
    inputPara >> invStr >> param.x1;
    inputPara >> invStr >> param.y0;
    inputPara >> invStr >> param.y1;
    inputPara >> invStr >> param.nx;
    inputPara >> invStr >> param.ny;
    inputPara >> invStr >> param.k;
    param.dimPk = (param.k + 1)*(param.k + 2)/2;
    inputPara >> invStr >> param.numGLP;
    inputPara >> invStr >> param.gamma;

    inputPara.close(); 
}