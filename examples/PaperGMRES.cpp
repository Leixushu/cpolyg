#include <cmath>
#include <iostream>
#include <iomanip>
#include "Meshes.h"
#include "Advection.h"
#include "Jacobian.h"
#include "TimeIntegration.h"
#include "Preconditioners.h"

using namespace std;

struct ExactGaussian : FnFunctor
{
    double theta;
    
    double operator()(double x, double y) const
    {
        double x0, y0;
        
        x0 = 0.5 - 0.15*cos(2*theta);
        y0 = 0.5 + 0.3*cos(theta)*sin(theta);
        
        return exp(-150*((x-x0)*(x-x0) + (y-y0)*(y-y0)));
    }
    
    ExactGaussian(double a_theta) : theta(a_theta) {};
};

double planarWave(double x, double y)
{
    int nx = 1;
    int ny = 1;
    return sin(nx*x*2*M_PI + ny*y*2*M_PI);
}

double zero(double x, double y)
{
    return 0;
}

void solveit(int meshType, double cfl, int deg)
{
    double h = 0.05;
    
    cout << "Creating mesh with h = " << h << endl;
    
    //h = 0.4;
    
    PolyMesh msh;
    switch(meshType)
    {
        case 0:
            cout << "(Naturally ordered) hexagons" << endl;
            msh = naturalOrderedHexRotated(h/sqrt(6.0), 1, 1);
            break;
        case 1:
            cout << "(Naturally ordered) squares" << endl;
            msh = naturalOrderedQuad(h*0.5*sqrt(sqrt(3)), 1, 1);
            break;
        case 2:
            cout << "(Naturally ordered) right triangles" << endl;
            msh = naturalOrderedTri(h*sqrt(sqrt(3)/2), 1, 1);
            break;
        case 3:
            cout << "(Naturally ordered) equilateral triangles" << endl;
            msh = naturalOrderedHoneycomb(h, 1, 1);
            break;
        case 4:
            cout << "Peturbed polygons" << endl;
            msh = perturbedQuadRectangle(h, 0.5, 1, 1);
            break;
        case 5:
            cout << "Peturbed triangles" << endl;
            msh = perturbedTriRectangle(h, 0.5, 1, 1);
            break;
    }
    
    cout << "Total degrees of freedom: " << msh.np*0.5*(deg+1)*(deg+2) << endl;
    
    msh.gnuplot();
    
    VecFnCallbackFunctor zeroFunctor(zero);
    BoundaryConditions bc = BoundaryConditions::dirichletConditions(msh, 
        &zeroFunctor);
    
    Advection eqn(msh, bc);
    
    MassMatrix M(msh, deg, eqn.nc);
    M.spy("plt/M.gnu");
    
    ExactGaussian exact(0);
    
    MeshFn f = MeshFn(msh, deg, 1);
    f.interp(exact);
    
    cout << "Computing Jacobian matrix" << endl;
    Jacobian B = eqn.jacobian(f, 0);
    
    MeshFn unp1 = f;
    
    int K = 3;
    int i;
    double dt = cfl*h/sqrt(2);
    
    cout << "Using dt = " << dt << endl;
    cout << "Computing total of " << K << " timesteps." << endl;
    
    B *= -dt;
    B += M;
    BlockILU0 pc(B);
    
    for (i = 0; i < K; i++)
    {
        if (i%10 == 0)
            cout << "Beginning timestep " << i << ", t=" << i*dt << endl;
        
        unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
        unp1 = B.solve(M.matvec(unp1), pc, kGMRESSolver);
    }
    unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
    
    cout << setprecision(20) << "Computed until final time t=" << i*dt << endl;
    
//     double l2err;
//     
//     exact.theta = i*dt;
//     l2err = unp1.L2Error(exact);
//     cout << "L^2 error = " << l2err << endl;
}

void doMeshes(double cfl, int deg)
{
    int i;
    
    cout << endl << endl << "CFL = " << cfl << endl << endl;
    
    solveit(3, cfl, deg);
    for (i = 0; i < 4; i++)
    {
        //solveit(i, cfl, deg);
    }
}

int main()
{
    int deg = 2;
    
    doMeshes(1, deg);
    doMeshes(2, deg);
    doMeshes(4, deg);
    
    return 0;
}
