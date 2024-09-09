#pragma once

# include "mesh.hpp"

class Output
{
    Mesh&        mesh_;

public:

    Output(Mesh& mesh);

    ~Output(){};

    void writeVtk(const std::string& filename, const array_3& data);


};