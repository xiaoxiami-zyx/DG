#pragma once

#include "mesh.hpp"
#include "field.hpp"
#include "base_function.hpp"
#include "output.hpp"

class Solver
{
    Mesh&           mesh_;
    Parameter&      param_;

    FieldCreator*    fieldCreator_ = nullptr;  // field creator
    BaseFunction*    baseFunction_ = nullptr;  // base function
    Gauss_Lobatto*   gl_ = nullptr;    // Gauss-Lobatto points

    Output*          output_ = nullptr;  // output

    std::map<std::string, Field> fields_;  // field container


public:
    Solver(Mesh& mesh, Parameter& param);

    ~Solver();

    void init();

    void createFields();

    void setInitialCondition();

};