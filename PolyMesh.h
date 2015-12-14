#pragma once

#include <vector>
#include <array>

struct PolyMesh
{
    int np;
    std::vector<std::vector<int>> p;
    std::vector<std::array<double, 2>> v;
    std::vector<std::array<double, 4>> bb;
    
    PolyMesh() {};
    PolyMesh(std::vector<std::array<double, 2>> pts);
};