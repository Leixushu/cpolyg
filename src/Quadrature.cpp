#include "Quadrature.h"

using namespace arma;

double Quadrature::polygonIntegral(const PolyMesh &msh, FnFunctor &cb, int pi, int deg)
{
    VecFnFunctor vecFunctor(cb);
    vec::fixed<1> integ = polygonIntegral(msh, vecFunctor, pi, deg);
    
    return integ(0);
}

mat Quadrature::polygonIntegral(const PolyMesh &msh, VecFunctor &cb, int pi, int deg)
{
    unsigned int i, k;
    double area, x1, x2, x3, y1, y2, y3, x, y;
    
    mat *qr;
    
    if (deg <= 2)
        qr = &tri2;
    else if (deg <= 4)
        qr = &tri4;
    else if (deg <= 6)
        qr = &tri6;
    else
        qr = &tri10;
    
    mat z = zeros<mat>(cb.n_rows, cb.n_cols);
    mat integ = zeros<mat>(cb.n_rows, cb.n_cols);
    
    // loop over triangles in triangulation of polygon
    for (i = 0; i < msh.tri[pi].triangles.size(); i++)
    {
        x1 = msh.tri[pi].p[i][0];
        y1 = msh.tri[pi].p[i][1];
        x2 = msh.tri[pi].p[i][2];
        y2 = msh.tri[pi].p[i][3];
        x3 = msh.tri[pi].p[i][4];
        y3 = msh.tri[pi].p[i][5];
        
        area = 0.5*fabs(x1*y2 + x2*y3 + x3*y1 - x1*y3 - x3*y2 - x2*y1);
        
        z *= 0;
        for (k = 0; k < qr->n_rows; k++)
        {
            x = x1*(1-(*qr)(k,0)-(*qr)(k,1))+x2*(*qr)(k,0)+x3*(*qr)(k,1);
            y = y1*(1-(*qr)(k,0)-(*qr)(k,1))+y2*(*qr)(k,0)+y3*(*qr)(k,1);
            
            z += cb(x, y)*(*qr)(k,2);
        }
        
        integ += z*area;
    }
    
    return integ;
}

double Quadrature::lineIntegral(const PolyMesh &msh, FnFunctor &cb, int a, int b, int deg)
{
    VecFnFunctor vecFunctor(cb);
    vec::fixed<1> integ = lineIntegral(msh, vecFunctor, a, b, deg);
    
    return integ(0);
}

mat Quadrature::lineIntegral(const PolyMesh &msh, VecFunctor &cb, int a, int b, int deg)
{
    unsigned int k;
    double x1, y1, x2, y2, x, y, length;
    
    mat *qr;

    if (deg <= 2)
        qr = &lin2;
    else if (deg <= 4)
        qr = &lin4;
    else if (deg <= 6)
        qr = &lin6;
    else
        qr = &lin10;
    
    mat integ = zeros<mat>(cb.n_rows, cb.n_cols);
    
    x1 = msh.v[a][0];
    y1 = msh.v[a][1];
    x2 = msh.v[b][0];
    y2 = msh.v[b][1];
    
    length = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
    
    for (k = 0; k < qr->n_rows; k++)
    {
        x = x1 + 0.5*((*qr)(k,0) + 1)*(x2 - x1);
        y = y1 + 0.5*((*qr)(k,0) + 1)*(y2 - y1);
        
        integ += (*qr)(k,1)*cb(x, y);
    }
    
    return integ*0.5*length;
}

mat Quadrature::lin2 = {{-0.5773502691896257, 1.0},
                         { 0.5773502691896257, 1.0}};
mat Quadrature::lin4 = {{-0.3399810435848563, 0.6521451548625461},
                         { 0.3399810435848563, 0.6521451548625461},
                         {-0.8611363115940526, 0.3478548451374538},
                         { 0.8611363115940526, 0.3478548451374538}};
mat Quadrature::lin6 = {{ 0.6612093864662645, 0.3607615730481386},
                         {-0.6612093864662645, 0.3607615730481386},
                         {-0.2386191860831969, 0.4679139345726910},
                         { 0.2386191860831969, 0.4679139345726910},
                         {-0.9324695142031521, 0.1713244923791704},
                         { 0.9324695142031521, 0.1713244923791704}};
mat Quadrature::lin10 = {{ 0.1488743389816312, 0.2955242247147529},
                          {-0.1488743389816312, 0.2955242247147529},
                          {-0.4333953941292472, 0.2692667193099963},
                          { 0.4333953941292472, 0.2692667193099963},
                          {-0.6794095682990244, 0.2190863625159820},
                          { 0.6794095682990244, 0.2190863625159820},
                          {-0.8650633666889845, 0.1494513491505806},
                          { 0.8650633666889845, 0.1494513491505806},
                          {-0.9739065285171717, 0.0666713443086881},
                          { 0.9739065285171717, 0.0666713443086881}};

mat Quadrature::tri2 = {{0.16666666666667, 0.16666666666667, 0.33333333333333},
                         {0.16666666666667, 0.66666666666667, 0.33333333333333},
                         {0.66666666666667, 0.16666666666667, 0.33333333333333}};
mat Quadrature::tri4 = {{0.44594849091597, 0.44594849091597, 0.22338158967801},
                         {0.44594849091597, 0.10810301816807, 0.22338158967801},
                         {0.10810301816807, 0.44594849091597, 0.22338158967801},
                         {0.09157621350977, 0.09157621350977, 0.10995174365532},
                         {0.09157621350977, 0.81684757298046, 0.10995174365532},
                         {0.81684757298046, 0.09157621350977, 0.10995174365532}};
mat Quadrature::tri6 = {{0.24928674517091, 0.24928674517091, 0.11678627572638},
                         {0.24928674517091, 0.50142650965818, 0.11678627572638},
                         {0.50142650965818, 0.24928674517091, 0.11678627572638},
                         {0.06308901449150, 0.06308901449150, 0.05084490637021},
                         {0.06308901449150, 0.87382197101700, 0.05084490637021},
                         {0.87382197101700, 0.06308901449150, 0.05084490637021},
                         {0.31035245103378, 0.63650249912140, 0.08285107561837},
                         {0.63650249912140, 0.05314504984482, 0.08285107561837},
                         {0.05314504984482, 0.31035245103378, 0.08285107561837},
                         {0.63650249912140, 0.31035245103378, 0.08285107561837},
                         {0.31035245103378, 0.05314504984482, 0.08285107561837},
                         {0.05314504984482, 0.63650249912140, 0.08285107561837}};
mat Quadrature::tri10 = {{0.33333333333333, 0.33333333333333, 0.09081799038275},
                          {0.48557763338366, 0.48557763338366, 0.03672595775647},
                          {0.48557763338366, 0.02884473323269, 0.03672595775647},
                          {0.02884473323269, 0.48557763338366, 0.03672595775647},
                          {0.10948157548504, 0.10948157548504, 0.04532105943553},
                          {0.10948157548504, 0.78103684902993, 0.04532105943553},
                          {0.78103684902993, 0.10948157548504, 0.04532105943553},
                          {0.30793983876412, 0.55035294182100, 0.07275791684542},
                          {0.55035294182100, 0.14170721941488, 0.07275791684542},
                          {0.14170721941488, 0.30793983876412, 0.07275791684542},
                          {0.55035294182100, 0.30793983876412, 0.07275791684542},
                          {0.30793983876412, 0.14170721941488, 0.07275791684542},
                          {0.14170721941488, 0.55035294182100, 0.07275791684542},
                          {0.24667256063990, 0.72832390459741, 0.02832724253106},
                          {0.72832390459741, 0.02500353476269, 0.02832724253106},
                          {0.02500353476269, 0.24667256063990, 0.02832724253106},
                          {0.72832390459741, 0.24667256063990, 0.02832724253106},
                          {0.24667256063990, 0.02500353476269, 0.02832724253106},
                          {0.02500353476269, 0.72832390459741, 0.02832724253106},
                          {0.06680325101220, 0.92365593358750, 0.00942166696373},
                          {0.92365593358750, 0.00954081540030, 0.00942166696373},
                          {0.00954081540030, 0.06680325101220, 0.00942166696373},
                          {0.92365593358750, 0.06680325101220, 0.00942166696373},
                          {0.06680325101220, 0.00954081540030, 0.00942166696373},
                          {0.00954081540030, 0.92365593358750, 0.00942166696373}};
