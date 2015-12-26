#include "Jacobian.h"
#include <cassert>

using namespace arma;

Jacobian::Jacobian(PolyMesh &a_msh, int a_deg, int a_nc) : msh(a_msh)
{
    int i;
    int j;
    int numBlocks;
    
    n_rows = msh.np;
    nc = a_nc;
    deg = a_deg;
    bl = (deg + 1)*(deg + 2)/2;
    
    mat zeroMat = zeros<mat>(bl, bl);
    
    numBlocks = 0;
    
    for (i = 0; i < n_rows; i++)
    {
        blocks.push_back(zeroMat);
        colIndices.push_back(i);
        rowBlock.push_back(numBlocks);
        numBlocks++;
        
        for (j = 0; j < msh.p2p[i].size(); j++)
        {
            blocks.push_back(zeroMat);
            colIndices.push_back(msh.p2p[i][j]);
            numBlocks++;
        }
    }
    
    rowBlock.push_back(numBlocks);
    nb = numBlocks;
}

MeshFn Jacobian::dot(const MeshFn &x)
{
    // right now multiple components not really implemented...!
    assert(nc == 1);
    
    MeshFn result(msh, deg, nc);
    mat component = x.a.tube(0, 0, n_rows-1, 0);
    vec b = matvec(vectorise(component, 1).t());
    
    result.a.tube(0, 0, n_rows-1, 0) = reshape(b, bl, n_rows).t();
    
    return result;
}

MeshFn Jacobian::solve(const MeshFn &b, Preconditioner &pc, Solver s)
{
    // right now multiple components not really implemented...!
    assert(nc == 1);
    
    MeshFn result(msh, deg, nc);
    mat component = b.a.tube(0, 0, n_rows-1, 0);
    vec bVec = vectorise(component, 1).t();
    vec x = zeros<vec>(bVec.n_rows);
    
    double tol = 1.e-12;
    int maxIt = 200;
    
    switch(s)
    {
        case kGMRESSolver:
            gmres(bVec, x, 20, tol, maxIt, pc);
            cout << "GMRES iterations: " << maxIt << endl;
            break;
        case kJacobiSolver:
            x = jacobi(bVec, tol, maxIt, pc);
            cout << "Jacobi iterations: " << maxIt << endl;
            break;
    }
    
    result.a.tube(0, 0, n_rows-1, 0) = reshape(x, bl, n_rows).t();
    
    return result;
}

Jacobian& Jacobian::operator +=(MassMatrix &M)
{
    int i;
    
    for (i = 0; i < n_rows; i++)
    {
        blocks[rowBlock[i]] += M.blocks[i];
    }
    
    return *this;
}
