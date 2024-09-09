#include "field.hpp"

Field::Field(size_t nx, size_t ny, size_t dimPk, size_t numGLP)
: nx_(nx), ny_(ny), dimPk_(dimPk), numGLP_(numGLP)
{
    data_.resize(boost::extents[nx_][ny_][dimPk_]);
    phi_.resize(boost::extents[nx_][ny_][numGLP_]);
    // std::fill(data_.data(), data_.data() + data_.num_elements(), 0.0);
    // std::fill(phi_.data(), phi_.data() + phi_.num_elements(), 0.0);
}

Field& Field::operator=(Field&& other) noexcept
{
    if (this != &other)
    {
        nx_ = other.nx_;
        ny_ = other.ny_;
        dimPk_ = other.dimPk_;
        numGLP_ = other.numGLP_;
        data_.resize(boost::extents[nx_][ny_][dimPk_]);
        phi_.resize(boost::extents[nx_][ny_][numGLP_]);
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
    // data_ = std::move(other.data_);  // 不支持move
    // phi_ = std::move(other.phi_);
    data_.resize(boost::extents[nx_][ny_][dimPk_]);
    phi_.resize(boost::extents[nx_][ny_][numGLP_]);
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
