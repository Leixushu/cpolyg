#include <cmath>
#include <iostream>
#include <iomanip>
#include "Meshes.h"
#include "EulerVortex.h"
#include "Jacobian.h"
#include "TimeIntegration.h"
#include "Preconditioners.h"

#define kGamma 1.4

int main(int argc, char ** argv)
{
    using namespace std;
    using namespace arma;
    
    int i;
    int deg = 1;
    double h = 0.5;
    
    cout << "Using h = " << h << endl;
    
    PolyMesh msh = hexRectangle(h, 20, 15);
    msh.gnuplot();
    
    MassMatrix M(msh, deg);
    M.spy("plt/M.gnu");
    
    EulerVortex eqn(msh, kGamma);
    
    MeshFn f = eqn.initialConditions(deg);
    MeshFn unp1 = f;
    
    f.gnuplot("plt/f.gnu");
    
//     int polyIdx = 11;
//     double h2 = 1.e-7;
//     Jacobian B = eqn.jacobian(f, 0);
//     B.spy("plt/B.gnu");
//     f.a(0,0,polyIdx) -= h2;
//     MeshFn b1 = eqn.assemble(f, 0);
//     f.a(0,0,polyIdx) += 2*h2;
//     MeshFn b2 = eqn.assemble(f, 0);
//     
//     cout.precision(11);
//     
//     ((vectorise(b2.a.slice(polyIdx)) - vectorise(b1.a.slice(polyIdx)))/(2*h2)).raw_print();
//     
//     cout << endl;
//     
//     B.blocks[B.rowBlock[polyIdx]].raw_print();
//     exit(0);
    
    BackwardEuler ti(M, eqn);
    
    double dt = h/10;
    
    int K = 10;
    
    cout << "Computing total of " << K << " timesteps." << endl;
    
    for (i = 0; i < K; i++)
    {
        unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
        cout << "Beginning timestep " << i + 1 << endl;
        unp1 = ti.advance(unp1, dt, i*dt);
    }
    unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
    
    eqn.exact->t = i*dt;
    vec error = unp1.L2Error(*eqn.exact);
    cout << setprecision(20) << "L^2 errors: " << endl;
    cout << "Density:    " << error(0) << endl;
    cout << "Velocity u: " << error(1) << endl;
    cout << "Velocity v: " << error(2) << endl;
    cout << "Energy:     " << error(3) << endl;
    
    return 0;
}
