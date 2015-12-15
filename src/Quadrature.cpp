#include "Quadrature.h"

using namespace std;

array<array<double, 2>, 2> lin2quadrature()
{
    array<array<double, 2>, 2> qr = 
    {{{{-0.5773502691896257, 1.0}},
      {{ 0.5773502691896257, 1.0}}}};
    
    return qr;
}

array<array<double, 2>, 4> lin4quadrature()
{
    array<array<double, 2>, 4> qr = 
    {{{{-0.3399810435848563, 0.6521451548625461}},
      {{ 0.3399810435848563, 0.6521451548625461}},
      {{-0.8611363115940526, 0.3478548451374538}},
      {{ 0.8611363115940526, 0.3478548451374538}}}};
    
    return qr;
}

array<array<double, 3>, 3> tri2quadrature()
{
    array<array<double, 3>, 3> qr = 
    {{{{0.16666666666667, 0.16666666666667, 0.33333333333333}},
      {{0.16666666666667, 0.66666666666667, 0.33333333333333}},
      {{0.66666666666667, 0.16666666666667, 0.33333333333333}}}};
    
    return qr;
}

array<array<double, 3>, 6> tri4quadrature()
{
    array<array<double, 3>, 6> qr = 
    {{{{0.44594849091597, 0.44594849091597, 0.22338158967801}},
      {{0.44594849091597, 0.10810301816807, 0.22338158967801}},
      {{0.10810301816807, 0.44594849091597, 0.22338158967801}},
      {{0.09157621350977, 0.09157621350977, 0.10995174365532}},
      {{0.09157621350977, 0.81684757298046, 0.10995174365532}},
      {{0.81684757298046, 0.09157621350977, 0.10995174365532}}}};
    
    return qr;
}

std::array<std::array<double, 2>, 2> Quadratures::lin2 = lin2quadrature();
std::array<std::array<double, 2>, 4> Quadratures::lin4 = lin4quadrature();

std::array<std::array<double, 3>, 3> Quadratures::tri2 = tri2quadrature();
std::array<std::array<double, 3>, 6> Quadratures::tri4 = tri4quadrature();
