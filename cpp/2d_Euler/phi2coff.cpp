#include "phi2coff.hpp"


Phi2Coff::Phi2Coff(const std::vector<scalar>& weight, const array_3& phiGauss, const std::vector<scalar>& mass)
:weight_(weight),phiGauss_(phiGauss),mm_(mass)
{
}

void Phi2Coff::computeCoff(Field& a)
{
    array_3& coff = a.data();
    array_4& phi = a.phi();

    auto [nx, ny, dimPk] = a.dataShape();
    auto [nx1, ny1, numGLP, numGLP2] = a.phiShape();

    assert(nx == nx1 && ny == ny1 && numGLP == numGLP2);

    array_3 temp(numGLP,numGLP,dimPk);
    for (label i1 = 0; i1 < numGLP; i1++)
    {
        for (label j1 = 0; j1 < numGLP; j1++)
        {
            scalar weight_ij = weight_[i1] * weight_[j1];
            for (label d = 0; d < dimPk; d++)
            {
                temp(i1,j1,d) = weight_ij * phiGauss_(i1,j1,d);
            }
        }
    }

    for(label i = 0; i < nx; i++){
    for(label j = 0; j < ny; j++){
    for(label d = 0; d < dimPk; d++){
        coff(i,j,d) = 0.0;
        for(label i1 = 0; i1 < numGLP; i1++){
        for(label j1 = 0; j1 < numGLP; j1++){
            coff(i,j,d) += temp(i1,j1,d) * phi(i,j,i1,j1);
        }
        }
        coff(i,j,d) *= 0.25 / mm_[d];
    }
    }
    }
}
