#pragma once
#include "param.hpp"
# include "mesh.hpp"
# include "gauss_Lobatto.hpp"
# include <boost/multi_array.hpp>

using array_3 = boost::multi_array<scalar, 3>;
using array_2 = boost::multi_array<scalar, 2>;

class field
{
    Mesh&        mesh_;
    Parameter&   param_;
    Gauss_Lobatto& glp_;  // Gauss-Lobatto 积分点和权重

    array_3      data_;   // 2d 系数场  nx ny dimPk 
    array_3      phi_;    // 2d 物理场  nx ny numGLP

    label        dimPk_;   // 系数场维数
    label        numGLP_;  // Gauss-Lobatto 积分点数
    label        nx_, ny_; // 网格数

public:
    field(Mesh& mesh, Parameter& param, Gauss_Lobatto& glp);

    ~field(){};

    const array_3& data() const { return data_; }

    const array_3& phi() const { return phi_; }
};