#include "field.hpp"

field::field(Mesh& mesh, Parameter& param, Gauss_Lobatto& glp)
    : mesh_(mesh), param_(param), glp_(glp), dimPk_(param_.dimPk),
    numGLP_(param_.numGLP), nx_(mesh_.nx()), ny_(mesh_.ny())
{
    data_.resize(boost::extents[nx_][ny_][dimPk_]);
    phi_.resize(boost::extents[nx_][ny_][numGLP_]);
}
