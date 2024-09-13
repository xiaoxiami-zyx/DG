#include "field.hpp"

Field::Field(label nx, label ny, label dimPk, label numGLP)
: nx_(nx), ny_(ny), dimPk_(dimPk), numGLP_(numGLP),
data_(nx_,ny_,dimPk_), phi_(nx_,ny_,numGLP_,numGLP_)
{
}

Field& Field::operator=(Field&& other) noexcept
{
    if (this != &other)
    {
        nx_ = other.nx_;
        ny_ = other.ny_;
        dimPk_ = other.dimPk_;
        numGLP_ = other.numGLP_;
        std::swap(data_, other.data_);
        std::swap(phi_, other.phi_);
    }
    return *this;
}

Field::Field(Field&& other) noexcept
{
    nx_ = other.nx_;
    ny_ = other.ny_;
    dimPk_ = other.dimPk_;
    numGLP_ = other.numGLP_;
    std::swap(data_, other.data_);
    std::swap(phi_, other.phi_);
}



FieldCreator::FieldCreator(Mesh& mesh, Parameter& param)
: mesh_(mesh), param_(param)
{
}

FieldCreator::~FieldCreator()
{
}

Field FieldCreator::createField() const
{
    return Field(mesh_.nx(), mesh_.ny(), param_.dimPk, param_.numGLP);
}
