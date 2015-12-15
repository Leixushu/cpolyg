#pragma once

#include <vector>
#include <array>
#include <cmath>

#include "Triangulation.h"
#include "Quadrature.h"

typedef double (*FnCallback)(double x, double y);

struct FnFunctor
{   
    FnCallback cb;
    
    FnFunctor() {};
    FnFunctor(FnCallback fn) : cb(fn) {};
    
    virtual double operator()(double x, double y)
    {
        return cb(x, y);
    }
};

struct PolyMesh
{
    int np;
    std::vector<std::vector<int>> p;
    std::vector<std::array<double, 2>> v;
    std::vector<std::array<double, 4>> bb;
    std::vector<std::vector<int>> p2p;
    std::vector<Triangulation> tri;
    
    PolyMesh() {};
    PolyMesh(std::vector<std::array<double, 2>> pts);
    
    int addVertex(std::array<double, 2> vertex);
    void computep2p();
    void computebb();
    void computeTriangulation();
    
    void getOutwardNormal(int i, int a, int b, double &x, double &y);
    
    inline void getLocalCoordinates(int p, double x_in, double y_in,
                                           double &x_out, double &y_out)
    {
        x_out = 2.0*(x_in - bb[p][0])/bb[p][2] - 1;
        y_out = 2.0*(y_in - bb[p][1])/bb[p][3] - 1;
    }
    
    inline double polygonIntegral(FnFunctor &cb, int p)
    {
        int i, k;
        double area, integ, x1, x2, x3, y1, y2, y3, x, y, z;
        auto &qr = Quadratures::tri4;
        
        integ = 0;
        
        // loop over triangles in triangulation of polygon
        for (i = 0; i < tri[p].triangles.size(); i++)
        {
            x1 = tri[p].p[i][0];
            y1 = tri[p].p[i][1];
            x2 = tri[p].p[i][2];
            y2 = tri[p].p[i][3];
            x3 = tri[p].p[i][4];
            y3 = tri[p].p[i][5];
            
            area = 0.5*fabs(x1*y2 + x2*y3 + x3*y1 - x1*y3 - x3*y2 - x2*y1);
            
            z = 0;
            for (k = 0; k < qr.size(); k++)
            {
                x = x1*(1-qr[k][0]-qr[k][1])+x2*qr[k][0]+x3*qr[k][1];
                y = y1*(1-qr[k][0]-qr[k][1])+y2*qr[k][0]+y3*qr[k][1];
                
                z += cb(x, y)*qr[k][2];
            }
            
            integ += z*area;
        }
        
        return integ;
    }
    
    inline double lineIntegral(FnFunctor &cb, int a, int b)
    {
        int k;
        double integ, x1, y1, x2, y2, x, y, length;
        auto &qr = Quadratures::lin4;
        
        integ = 0;
        
        x1 = v[a][0];
        y1 = v[a][1];
        x2 = v[b][0];
        y2 = v[b][1];
        
        length = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
        
        for (k = 0; k < qr.size(); k++)
        {
            x = x1 + 0.5*(qr[k][0] + 1)*(x2 - x1);
            y = y1 + 0.5*(qr[k][0] + 1)*(y2 - y1);
            
            integ += qr[k][1]*cb(x, y);
        }
        
        return integ*0.5*length;
    }
    
    double integrate(FnFunctor &cb);
};