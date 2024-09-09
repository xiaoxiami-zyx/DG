#pragma once
#include "param.hpp"
# include "mesh.hpp"


class FieldCreator; // 前置声明
class Field
{
    array_3      data_;   // 2d 系数场  nx ny dimPk  每个cell多项式系数 
    array_3      phi_;    // 2d 物理场  nx ny numGLP 每个cell物理场值， 一个cell内有多个积分点 

    size_t        nx_, ny_; // 网格数
    size_t        dimPk_;   // 多项式次数
    size_t        numGLP_;  // Gauss-Lobatto 积分点数

private:
    Field(size_t nx, size_t ny, size_t dimPk, size_t numGLP);  // 构造函数私有化

    friend class FieldCreator; // 友元类

public:
    ~Field(){};

    // 禁用拷贝构造函数
    Field(const Field&) = delete;

    // 禁用拷贝赋值运算符
    Field& operator=(const Field&) = delete;

    // 移动赋值运算符
    Field &operator=(Field &&other) noexcept;

    // 移动构造函数
    Field(Field &&other) noexcept;


    constexpr scalar& operator[](size_t i,size_t j,size_t k) noexcept { return data_[i][j][k];}

    constexpr const scalar& operator[](size_t i, size_t j, size_t k) const noexcept { return data_[i][j][k]; }

    array_3& data() noexcept { return data_; }

    array_3& phi() noexcept { return phi_; }

    std::tuple<size_t, size_t ,size_t> dataShape() const noexcept { return std::make_tuple(nx_, ny_, dimPk_); }

    std::tuple<size_t, size_t ,size_t> phiShape() const noexcept { return std::make_tuple(nx_, ny_, numGLP_); }

};

class FieldCreator
{
    Mesh&        mesh_;
    Parameter&   param_;

public:
    FieldCreator(Mesh& mesh, Parameter& param);

    ~FieldCreator();

    Field createField() const;

};