#include "base_function.hpp"

BaseFunction::BaseFunction(label dimPk, label numGLP, Gauss_Lobatto& glp, scalar dx, scalar dy)
: dimPk_(dimPk), numGLP_(numGLP),
phiGauss_(numGLP_, numGLP_, dimPk_), 
phiGaussLL_(numGLP_, numGLP_, dimPk_, 3),
phiGaussR_(numGLP_, dimPk_), 
phiGaussU_(numGLP_, dimPk_), 
phiGaussD_(numGLP_, dimPk_), 
phiGaussL_(numGLP_, dimPk_),
phiGauss_dx_(numGLP_, numGLP_, dimPk_),
phiGauss_dy_(numGLP_, numGLP_, dimPk_)
{
    phiLD_.resize(dimPk_);
    phiLU_.resize(dimPk_);
    phiRD_.resize(dimPk_);
    phiRU_.resize(dimPk_);
    mm_.resize(dimPk_);

    phiRU_[0] = basefunc_0(1.0, 1.0);
    phiRU_[1] = basefunc_1(1.0, 1.0);
    phiRU_[2] = basefunc_2(1.0, 1.0);
    phiRU_[3] = basefunc_3(1.0, 1.0);
    phiRU_[4] = basefunc_4(1.0, 1.0);
    phiRU_[5] = basefunc_5(1.0, 1.0);

    phiLU_[0] = basefunc_0(-1.0, 1.0);
    phiLU_[1] = basefunc_1(-1.0, 1.0);
    phiLU_[2] = basefunc_2(-1.0, 1.0);
    phiLU_[3] = basefunc_3(-1.0, 1.0);
    phiLU_[4] = basefunc_4(-1.0, 1.0);
    phiLU_[5] = basefunc_5(-1.0, 1.0);

    phiRD_[0] = basefunc_0(1.0, -1.0);
    phiRD_[1] = basefunc_1(1.0, -1.0);
    phiRD_[2] = basefunc_2(1.0, -1.0);
    phiRD_[3] = basefunc_3(1.0, -1.0);
    phiRD_[4] = basefunc_4(1.0, -1.0);

    phiLD_[0] = basefunc_0(-1.0, -1.0);
    phiLD_[1] = basefunc_1(-1.0, -1.0);
    phiLD_[2] = basefunc_2(-1.0, -1.0);
    phiLD_[3] = basefunc_3(-1.0, -1.0);
    phiLD_[4] = basefunc_4(-1.0, -1.0);
    phiLD_[5] = basefunc_5(-1.0, -1.0);

    mm_[0] = 1.0;      mm_[1] = 1.0/3.0;   mm_[2] = 1.0/3.0; 
    mm_[3] = 4.0/45.0; mm_[4] = 1.0/9.0;   mm_[5] = 4.0/45.0;

    const std::vector<scalar>& lambda = glp.lambda();
    const std::vector<scalar>& lambdaL = glp.lambdaL();

    for(label i = 0; i < numGLP_; i++)
    for(label j = 0; j < numGLP_; j++)
    {
        phiGauss_(i,j,0) = basefunc_0(lambda[i], lambda[j]);
        phiGauss_(i,j,1) = basefunc_1(lambda[i], lambda[j]);
        phiGauss_(i,j,2) = basefunc_2(lambda[i], lambda[j]);
        phiGauss_(i,j,3) = basefunc_3(lambda[i], lambda[j]);
        phiGauss_(i,j,4) = basefunc_4(lambda[i], lambda[j]);
        phiGauss_(i,j,5) = basefunc_5(lambda[i], lambda[j]);

        phiGaussLL_(i,j,0,0) = basefunc_0(lambdaL[i], lambda[j]);
        phiGaussLL_(i,j,0,1) = basefunc_0(lambda[i], lambdaL[j]);
        phiGaussLL_(i,j,0,2) = basefunc_0(lambda[i], lambda[j]);
        phiGaussLL_(i,j,1,0) = basefunc_1(lambdaL[i], lambda[j]);
        phiGaussLL_(i,j,1,1) = basefunc_1(lambda[i], lambdaL[j]);
        phiGaussLL_(i,j,1,2) = basefunc_1(lambda[i], lambda[j]);
        phiGaussLL_(i,j,2,0) = basefunc_2(lambdaL[i], lambda[j]);
        phiGaussLL_(i,j,2,1) = basefunc_2(lambda[i], lambdaL[j]);
        phiGaussLL_(i,j,2,2) = basefunc_2(lambda[i], lambda[j]);
        phiGaussLL_(i,j,3,0) = basefunc_3(lambdaL[i], lambda[j]);
        phiGaussLL_(i,j,3,1) = basefunc_3(lambda[i], lambdaL[j]);
        phiGaussLL_(i,j,3,2) = basefunc_3(lambda[i], lambda[j]);
        phiGaussLL_(i,j,4,0) = basefunc_4(lambdaL[i], lambda[j]);
        phiGaussLL_(i,j,4,1) = basefunc_4(lambda[i], lambdaL[j]);
        phiGaussLL_(i,j,4,2) = basefunc_4(lambda[i], lambda[j]);
        phiGaussLL_(i,j,5,0) = basefunc_5(lambdaL[i], lambda[j]);
        phiGaussLL_(i,j,5,1) = basefunc_5(lambda[i], lambdaL[j]);
        phiGaussLL_(i,j,5,2) = basefunc_5(lambda[i], lambda[j]);


        phiGauss_dx_(i,j,0) = basefunc_dx_0(lambda[i], lambda[j], dx / 2.0, dy / 2.0);
        phiGauss_dx_(i,j,1) = basefunc_dx_1(lambda[i], lambda[j], dx / 2.0, dy / 2.0);
        phiGauss_dx_(i,j,2) = basefunc_dx_2(lambda[i], lambda[j], dx / 2.0, dy / 2.0);
        phiGauss_dx_(i,j,3) = basefunc_dx_3(lambda[i], lambda[j], dx / 2.0, dy / 2.0);
        phiGauss_dx_(i,j,4) = basefunc_dx_4(lambda[i], lambda[j], dx / 2.0, dy / 2.0);
        phiGauss_dx_(i,j,5) = basefunc_dx_5(lambda[i], lambda[j], dx / 2.0, dy / 2.0);
        phiGauss_dy_(i,j,0) = basefunc_dy_0(lambda[i], lambda[j], dx / 2.0, dy / 2.0);
        phiGauss_dy_(i,j,1) = basefunc_dy_1(lambda[i], lambda[j], dx / 2.0, dy / 2.0);
        phiGauss_dy_(i,j,2) = basefunc_dy_2(lambda[i], lambda[j], dx / 2.0, dy / 2.0);
        phiGauss_dy_(i,j,3) = basefunc_dy_3(lambda[i], lambda[j], dx / 2.0, dy / 2.0);
        phiGauss_dy_(i,j,4) = basefunc_dy_4(lambda[i], lambda[j], dx / 2.0, dy / 2.0);
        phiGauss_dy_(i,j,5) = basefunc_dy_5(lambda[i], lambda[j], dx / 2.0, dy / 2.0);
    }

    for(label i = 0; i < numGLP_; i++)
    {
        phiGaussR_(i,0) = basefunc_0(1.0, lambda[i]);
        phiGaussR_(i,1) = basefunc_1(1.0, lambda[i]);
        phiGaussR_(i,2) = basefunc_2(1.0, lambda[i]);
        phiGaussR_(i,3) = basefunc_3(1.0, lambda[i]);
        phiGaussR_(i,4) = basefunc_4(1.0, lambda[i]);
        phiGaussR_(i,5) = basefunc_5(1.0, lambda[i]);

        phiGaussU_(i,0) = basefunc_0(lambda[i], 1.0);
        phiGaussU_(i,1) = basefunc_1(lambda[i], 1.0);
        phiGaussU_(i,2) = basefunc_2(lambda[i], 1.0);
        phiGaussU_(i,3) = basefunc_3(lambda[i], 1.0);
        phiGaussU_(i,4) = basefunc_4(lambda[i], 1.0);
        phiGaussU_(i,5) = basefunc_5(lambda[i], 1.0);

        phiGaussD_(i,0) = basefunc_0(lambda[i], -1.0);
        phiGaussD_(i,1) = basefunc_1(lambda[i], -1.0);
        phiGaussD_(i,2) = basefunc_2(lambda[i], -1.0);
        phiGaussD_(i,3) = basefunc_3(lambda[i], -1.0);
        phiGaussD_(i,4) = basefunc_4(lambda[i], -1.0);
        phiGaussD_(i,5) = basefunc_5(lambda[i], -1.0);

        phiGaussL_(i,0) = basefunc_0(-1.0, lambda[i]);
        phiGaussL_(i,1) = basefunc_1(-1.0, lambda[i]);
        phiGaussL_(i,2) = basefunc_2(-1.0, lambda[i]);
        phiGaussL_(i,3) = basefunc_3(-1.0, lambda[i]);
        phiGaussL_(i,4) = basefunc_4(-1.0, lambda[i]);
        phiGaussL_(i,5) = basefunc_5(-1.0, lambda[i]);
    }

}