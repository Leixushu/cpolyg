#pragma once
#include <map>
#include <memory>
#include <armadillo>
#include "Functors.h"
#include "PolyMesh.h"

enum BoundaryType
{
    kDirichletCondition,
    kPeriodicCondition
};

struct BoundaryConditions
{
    // Boundary conditions can be time dependent
    double t = 0;
    
    const PolyMesh &msh;
    
    std::map<int, BoundaryType> bcTypes;
    std::map<int, VecFunctor *> dirichlet;
    
    arma::mat boundaryValue(double x, double y, arma::mat U, int b);
    
    BoundaryConditions(PolyMesh &a_msh) : msh(a_msh) { }
    
    static BoundaryConditions dirichletConditions(PolyMesh &a_mhs, 
        VecFunctor *fn);
    static BoundaryConditions periodicConditions(PolyMesh &a_mhs);
};
