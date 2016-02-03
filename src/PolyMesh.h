#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <armadillo>
#include <cassert>

#include "Triangulation.h"
#include "Quadrature.h"
#include "Functors.h"

#define kEPS 1.e-10

struct PolyMesh
{
    int np;
    std::vector<std::vector<int> > p;
    std::vector<std::array<double, 2> > v;
    std::vector<std::array<double, 4> > bb;
    std::vector<std::vector<int> > p2p;
    std::vector<Triangulation> tri;
    
    PolyMesh() {};
    PolyMesh(std::vector<std::array<double, 2> > pts, double width = 1, double height = 1);
    static PolyMesh triangulate(std::vector<std::array<double, 2> > pts);
    
    int addVertex(std::array<double, 2> vertex);
    void computep2p();
    void computebb();
    void computeTriangulation();
    
    void gnuplot();
    
    void getOutwardNormal(int i, int a, int b, double &x, double &y);
    
    inline void getLocalCoordinates(int pi, double x_in, double y_in,
                                            double &x_out, double &y_out) const
    {
        x_out = 2.0*(x_in - bb[pi][0])/bb[pi][2] - 1;
        y_out = 2.0*(y_in - bb[pi][1])/bb[pi][3] - 1;
    }
    
    inline void getPeriodicCoordinates(int p1, int neighbor, double xIn, double yIn,
                                       double &xOut, double &yOut) const
    {
        int e1, e2, p2;
        double x1, y1, x2, y2, x3, y3, x4, y4, xNew, yNew;
        
        p2 = -neighbor - 1;
        for (e1 = 0; e1 < p2p[p1].size(); e1++)
        {
            if (p2p[p1][e1] == neighbor)
            {
                break;
            }
        }
        assert(e1 < p2p[p1].size());
        
        for (e2 = 0; e2 < p2p[p2].size(); e2++)
        {
            if (p2p[p2][e2] == -p1 - 1)
            {
                break;
            }
        }
        assert(e2 < p2p[p2].size());
        
        x1 = v[p[p1][e1]][0];
        y1 = v[p[p1][e1]][1];
        
        x2 = v[p[p2][e2]][0];
        y2 = v[p[p2][e2]][1];
        
        x3 = v[p[p1][e1+1]][0];
        y3 = v[p[p1][e1+1]][1];
        
        x4 = v[p[p2][e2+1]][0];
        y4 = v[p[p2][e2+1]][1];
        
        xNew = (x4*y1 - x2*y3)/(x3*y1 - x1*y3)*xIn + (x2*x3 - x1*x4)/(x3*y1 - x1*y3)*yIn;
        yNew = (y2*y3 - y1*y4)/(-(x3*y1) + x1*y3)*xIn + (x3*y2 - x1*y4)/(x3*y1 - x1*y3)*yIn;
        
        xOut = 2.0*(xNew - bb[p2][0])/bb[p2][2] - 1;
        yOut = 2.0*(yNew - bb[p2][1])/bb[p2][3] - 1;
    }
    
    inline double polygonIntegral(FnFunctor &cb, int pi)
    {
        unsigned int i, k;
        double area, integ, x1, x2, x3, y1, y2, y3, x, y, z;
        auto &qr = Quadratures::tri6;
        
        integ = 0;
        
        // loop over triangles in triangulation of polygon
        for (i = 0; i < tri[pi].triangles.size(); i++)
        {
            x1 = tri[pi].p[i][0];
            y1 = tri[pi].p[i][1];
            x2 = tri[pi].p[i][2];
            y2 = tri[pi].p[i][3];
            x3 = tri[pi].p[i][4];
            y3 = tri[pi].p[i][5];
            
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
    
    inline arma::mat polygonIntegral(VecFunctor &cb, int pi)
    {
        unsigned int i, k;
        double area, x1, x2, x3, y1, y2, y3, x, y;
        auto &qr = Quadratures::tri6;
        
        arma::mat z = arma::zeros<arma::mat>(cb.n_rows, cb.n_cols);
        arma::mat integ = arma::zeros<arma::mat>(cb.n_rows, cb.n_cols);
        
        // loop over triangles in triangulation of polygon
        for (i = 0; i < tri[pi].triangles.size(); i++)
        {
            x1 = tri[pi].p[i][0];
            y1 = tri[pi].p[i][1];
            x2 = tri[pi].p[i][2];
            y2 = tri[pi].p[i][3];
            x3 = tri[pi].p[i][4];
            y3 = tri[pi].p[i][5];
            
            area = 0.5*fabs(x1*y2 + x2*y3 + x3*y1 - x1*y3 - x3*y2 - x2*y1);
            
            z *= 0;
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
        unsigned int k;
        double integ, x1, y1, x2, y2, x, y, length;
        auto &qr = Quadratures::lin6;
        
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
    
    inline arma::mat lineIntegral(VecFunctor &cb, int a, int b)
    {
        unsigned int k;
        double x1, y1, x2, y2, x, y, length;
        auto &qr = Quadratures::lin6;
        
        arma::mat integ = arma::zeros<arma::mat>(cb.n_rows, cb.n_cols);
        
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
