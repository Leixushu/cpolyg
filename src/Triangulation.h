#pragma once

#include <vector>
#include <array>

extern "C" {
#include "triangle/triangle.h"
}

struct Triangulation
{
    std::vector<std::array<double, 2> > points;
    std::vector<std::array<int, 3> > triangles;
    std::vector<std::array<double, 6> > p;
    
    Triangulation() { }
    Triangulation(std::vector<std::array<double, 2> > pts);
    void doTriangulation();
};
