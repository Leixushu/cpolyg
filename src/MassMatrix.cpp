#include "MassMatrix.h"
#include "Legendre.h"

using namespace arma;
using namespace std;

struct ProductFunctor : FnFunctor
{
    vec phi, psi;
    int i;
    int m;
    PolyMesh &msh;
    
    ProductFunctor(PolyMesh &m) : msh(m) {};
    
    double operator()(double x, double y)
    {
        double xx, yy;
        msh.getLocalCoordinates(i, x, y, xx, yy);
        
        return Leg2D(xx, yy, m, phi)*Leg2D(xx, yy, m, psi);
    };
};

MassMatrix::MassMatrix(PolyMesh &m, int d) : msh(m), deg(d)
{
    int i, j, k;
    int basisSize = (deg+1)*(deg+2)/2;
    int N;
    double integ;
    ProductFunctor prod(msh);
    
    prod.phi = zeros<vec>(basisSize);
    prod.psi = zeros<vec>(basisSize);
    
    prod.m = deg+1;
    
    N = msh.np*basisSize;
    
    matrix  = sp_mat(N, N);
    
    for (i = 0; i < msh.np; i++)
    {
        prod.i = i;
        
        for (j = 0; j < basisSize; j++)
        {
            prod.phi[j] = 1.0;
            for (k = 0; k < basisSize; k++)
            {
                prod.psi[k] = 1.0;
                
                integ = msh.polygonIntegral(prod, i);
                
                matrix(i*basisSize + j, i*basisSize + k) = integ;
                matrix(i*basisSize + k, i*basisSize + j) = integ;
                
                prod.psi[k] = 0.0;
            }
            prod.phi[j] = 0.0;
        }
    }
}

MeshFn MassMatrix::solve(const MeshFn &fn) const
{
    int basisSize = (deg+1)*(deg+2)/2;
    MeshFn result(msh, deg, 1);
    
    mat firstComp = fn.a.tube(0, 0, result.a.n_rows - 1, 0);
    vec b = vectorise(firstComp, 1).t();
    vec x = arma::spsolve(matrix, b);
    
    result.a.tube(0, 0, result.a.n_rows - 1, 0) = reshape(x, basisSize, msh.np).t();
    
    return result;
}

void MassMatrix::spy(std::string filename)
{
    int i, j;
    ofstream file;
    file.open(filename);
    
    for (i = 0; i < matrix.n_rows; i++)
    {
        for (j = 0; j < matrix.n_cols; j++)
        {
            if (matrix(i, j) != 0.0)
            {
                file << i << "\t" << j << "\t" << matrix(i, j) << endl;
            }
        }
    }
    
    file.close();
}
