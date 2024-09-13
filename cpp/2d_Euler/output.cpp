#include "output.hpp"

Output::Output(Mesh& mesh)
    : mesh_(mesh)
{
}

void Output::writeVtk(const std::string& filename, const array_4& data)
{
    array_2 temp(mesh_.nx(), mesh_.ny());
    for (label i = 0; i < mesh_.nx(); i++) {
        for (label j = 0; j < mesh_.ny(); j++) {
            temp(i,j) = data(i,j,0,0);
        }
    }
    writeVtk(filename, temp);
}

void Output::writeVtk(const std::string& filename, const array_3& data)
{
    array_2 temp(mesh_.nx(), mesh_.ny());
    for (label i = 0; i < mesh_.nx(); i++) {
        for (label j = 0; j < mesh_.ny(); j++) {
            temp(i,j) = data(i,j,0);
        }
    }
    writeVtk(filename, temp);
}




void Output::writeVtk(const std::string& filename, const array_2& data)
{
    std::vector<scalar> x = mesh_.x();
    std::vector<scalar> y = mesh_.y();

    label nx = mesh_.nx() + 1;
    label ny = mesh_.ny() + 1;

    std::ofstream file(filename+".vtk");

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    file << "# vtk DataFile Version 3.0\n";
    file << "Structured Grid Example\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_GRID\n";
    file << "DIMENSIONS " << nx << " " << ny << " " << 1 << "\n";
    file << "POINTS " << nx * ny << " float\n";

    for (label i = 0; i < nx; ++i) {
        for (label j = 0; j < ny; ++j) {
            file << x[i] << " " << y[j] << " " << 0.0 << "\n";
        }
    }

    file << "CELL_DATA " << (nx - 1) * (ny - 1) * 1 << "\n";
    file << "SCALARS cell_data float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (label i = 0; i < nx - 1; ++i) {
        for (label j = 0; j < ny - 1; ++j) {
                file << data(i,j) << "\n";
        }
    }

}
