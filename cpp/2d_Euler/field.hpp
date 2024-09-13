#pragma once
#include "param.hpp"
# include "mesh.hpp"


class FieldCreator; // 前置声明
class Field
{
    label        nx_, ny_; // 网格数
    label        dimPk_;   // 多项式次数
    label        numGLP_;  // Gauss-Lobatto 积分点数

    array_3      data_;   // 2d 系数场  nx ny dimPk  每个cell多项式系数 
    array_4      phi_;    // 2d 物理场  nx ny numGLP numGLP 每个cell物理场值， 一个cell内有多个积分点 

private:
    Field(label nx, label ny, label dimPk, label numGLP);  // 构造函数私有化

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


    // constexpr scalar& operator[](label i,label j,label k) noexcept { return data_[i][j][k];}

    // constexpr const scalar& operator[](label i, label j, label k) const noexcept { return data_[i][j][k]; }

    array_3& data() noexcept { return data_; }

    array_4& phi() noexcept { return phi_; }

    std::tuple<label, label ,label> dataShape() const noexcept { return std::make_tuple(nx_, ny_, dimPk_); }

    std::tuple<label, label ,label, label> phiShape() const noexcept { return std::make_tuple(nx_, ny_, numGLP_,numGLP_); }

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