#pragma once

#include <array>

struct Quadratures
{
    static std::array<std::array<double, 2>, 2> lin2;
    static std::array<std::array<double, 2>, 4> lin4;
    
    static std::array<std::array<double, 3>, 3> tri2;
    static std::array<std::array<double, 3>, 6> tri4;
};
