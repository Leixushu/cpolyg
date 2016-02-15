#pragma once

#include <vector>
#include <array>
#include <map>
#include <armadillo>

#include "Triangulation.h"
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
    int findVertex(double x, double y);
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
    
    virtual ~PolyMesh() {};
};

struct PeriodicMesh : PolyMesh
{
    struct PeriodicBC
    {
        int p1, p2;
        int a1, b1, a2, b2;
    };
    
    std::vector<PeriodicBC> bc;
    
    void getPeriodicCoordinates(int p1, int p2, double x_in, double y_in,
                                double &x_out, double &y_out) const;
};
