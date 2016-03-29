#include <cmath>
#include "Meshes.h"
#include "Advection.h"
#include "TimeIntegration.h"
#include "Preconditioners.h"

// Initial conditions
double gaussian(double x, double y)
{
    return exp(-150*((x-0.25)*(x-0.25) + (y-0.5)*(y-0.5)));
}

double zero(double x, double y)
{
    return 0;
}

int main(int argc, char ** argv)
{
    using namespace std;
    // specify the degree of polynomials to use
    int deg = 4;
    // specify the size of each element in the mesh
    double h = 0.1;
    // specify the time step
    double dt = h*0.5;
    
    // create the mesh, in this case hexagonal tessellation
    // of the unit square
    PolyMesh msh = hexUnitSquare(h);
    
    // specify the advection equation
    // with zero Dirichlet boundary conditions
    VecFnCallbackFunctor zeroFunctor(zero);
    BoundaryConditions bc = 
        BoundaryConditions::dirichletConditions(msh, &zeroFunctor);
    Advection eqn(msh, bc);
    
    // compute the mass matrix
    MassMatrix M(msh, deg);
    
    MeshFn f = MeshFn(msh, gaussian, deg);
    
    int K = 2;
    
    DIRK3 tiDIRK(M, eqn);
    f = MeshFn(msh, gaussian, deg);
    for (int i = 0; i < K; i++)
    {
        cout << "Beginning timestep " << i << endl;
        f = tiDIRK.advance(f, dt, i*dt);
        f.gnuplot("plt/dirk" + to_string(i) + ".gnu");
    }
    
    // create the time integration object
    IRK3 ti(M, eqn);
    
    // create a mesh function by interpolating the 
    // initial conditions defined by the function above
    f = MeshFn(msh, gaussian, deg);
    for (int i = 0; i < K; i++)
    {
        cout << "Beginning timestep " << i << endl;
        f = ti.advance(f, dt, i*dt);
        f.gnuplot("plt/irk1" + std::to_string(i) + ".gnu");
    }
    
    f = MeshFn(msh, gaussian, deg);
    for (int i = 0; i < K; i++)
    {
        cout << "Beginning timestep " << i << endl;
        f = ti.newAdvance(f, dt, i*dt);
        f.gnuplot("plt/irk2" + to_string(i) + ".gnu");
    }
    
    return 0;
}
