#include <cmath>
#include <iostream>
#include <iomanip>
#include "Meshes.h"
#include "Advection.h"
#include "Jacobian.h"
#include "TimeIntegration.h"
#include "Preconditioners.h"

using namespace std;
using namespace arma;

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

double solveIt(int deg, double h)
{
    PolyMesh msh = quadUnitSquare(h);
    msh.gnuplot();
    
    MassMatrix M(msh, deg);
    
    Advection eqn(msh);
    ExactGaussian exact(0);
    
    MeshFn f = MeshFn(msh, deg, 1);
    f.interp(exact);
    
    MeshFn unp1 = f;
    
    int K = 10000;
    int i;
    double dt = M_PI/4/K;
    
    cout << "Using degree " << deg << ", h = " << h << ", dt = " << dt << endl;
    //cout << "Computing total of " << K << " timesteps." << endl;
    
    Jacobian B = eqn.jacobian(f, 0);
    
    for (i = 0; i < K; i++)
    {
        if (i%500 == 0)
        {
            cout << ".";
        }
        
        MeshFn k1 = dt*M.solve(B.dot(unp1));
        MeshFn k2 = dt*M.solve(B.dot(unp1 + 0.5*k1));
        MeshFn k3 = dt*M.solve(B.dot(unp1 + 0.5*k2));
        MeshFn k4 = dt*M.solve(B.dot(unp1 + k3));
        
        unp1 += (k1 + 2*k2 + 2*k3 + k4)/6;
    }
    cout << endl;
    
    double l2err;
    exact.theta = i*dt;
    
    MeshFn theExact = MeshFn(msh, deg, 1);
    theExact.interp(exact);
    theExact.gnuplot("plt/exact.gnu");
    
    l2err = unp1.L2Error(exact);
    cout << setprecision(20) << "L^2 error = " << l2err << endl << endl;
    
    return l2err;
}

void convergenceTest(int deg)
{
    double ratio = 0.75;
    int numTests = 4;
    vec e(numTests);
    int i;
    
    cout << "### Degree " << deg << " ###########################################"<<endl;
    for (i = 0; i < numTests; i++)
    {
        e(i) = solveIt(deg, 0.1*pow(ratio, i));
    }
    
    cout << "Orders:\t\t";
    for (i = 1; i < numTests; i++)
    {
        cout << log(e(i-1)/e(i))/(-log(ratio)) << "\t\t";
    }
    cout << endl;
    
    e.save("deg" + to_string(deg) + ".dat", csv_ascii);   
}

int main(int argc, char ** argv)
{
//     convergenceTest(0);
//     convergenceTest(1);
//     convergenceTest(2);
//     convergenceTest(3);
//     convergenceTest(4);
//     convergenceTest(5);
    convergenceTest(6);
    
    return 0;
}
