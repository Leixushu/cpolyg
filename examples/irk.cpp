#include "Advection.h"
#include "Meshes.h"
#include "Preconditioners.h"
#include "TimeIntegration.h"
#include "EulerSetups.h"
#include <cmath>
#include <iomanip>

#define EULER_EQN

const double kGamma = 1.4;

// Initial conditions
double gaussian(double x, double y) {
    return exp(-150 * ((x - 0.25) * (x - 0.25) + (y - 0.5) * (y - 0.5)));
}

double zero(double x, double y) { return 0; }

int main(int argc, char **argv) {
    using namespace std;
    
    #ifdef EULER_EQN
        int deg = 4;
        double h = 3;
        double dt = 0.075;
        PolyMesh msh = quadRectangle(h, 20, 15);
    #else
        int deg = 4;
        double h = 0.1;
        double dt = h * 0.5;
        PolyMesh msh = hexUnitSquare(h);
    #endif
    
    msh.gnuplot();

    string DOFStr = "%    Total DOF:  " +
                    to_string(int(msh.np * (deg + 1) * (deg + 2) * 0.5));
    string blockStr =
        "%    Block Size: " + to_string(int((deg + 1) * (deg + 2) / 2));

    cout << endl
         << "%--------------------------------------------------%" << endl
         << DOFStr << setw(52 - DOFStr.length()) << "%" << endl
         << blockStr << setw(52 - blockStr.length()) << "%" << endl
         << "%--------------------------------------------------%" << endl
         << endl;
    
    #ifdef EULER_EQN
        VortexSolution exactSolution(kGamma);
        Euler eqn(msh, 
            BoundaryConditions::dirichletConditions(msh, &exactSolution), kGamma);
        MeshFn f = exactSolution.interpolated(msh, deg);;
    #else
        VecFnCallbackFunctor zeroFunctor(zero);
        BoundaryConditions bc =
        BoundaryConditions::dirichletConditions(msh, &zeroFunctor);
        Advection eqn(msh, bc);
        MeshFn f = MeshFn(msh, gaussian, deg);
    #endif
    
    MeshFn u;

    // compute the mass matrix
    MassMatrix M(msh, deg, eqn.nc);

    int K = 1;

    // create the time integration object
    IRK ti(2, M, eqn);
    
    u = f;
    for (int i = 0; i < K; i++) {
        cout << "Beginning timestep " << i << endl;
        u = ti.newAdvance(u, dt, i * dt);
        u.gnuplot("plt/irk2" + to_string(i) + ".gnu");
    }
    
    u = f;
    for (int i = 0; i < K; i++) {
        cout << "Beginning timestep " << i << endl;
        u = ti.advance(u, dt, i * dt);
        u.gnuplot("plt/irk1" + std::to_string(i) + ".gnu");
    }

    DIRK3 tiDIRK(M, eqn);
    u = f;
    for (int i = 0; i < K; i++) {
        cout << "Beginning timestep " << i << endl;
        u = tiDIRK.advance(u, dt, i * dt);
        u.gnuplot("plt/dirk" + to_string(i) + ".gnu");
    }

    return 0;
}
