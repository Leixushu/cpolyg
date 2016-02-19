#include "Equation.h"
#include "Quadrature.h"

using namespace arma;
using namespace std;

typedef arma::mat (Equation::*Integrand)(double x, double y);
    
struct Equation::IntegrandFunctor : VecFunctor
{
    Equation &eqn;
    Integrand integ;
    
    IntegrandFunctor(Equation &a_eqn, Integrand a_integ, int a_n_rows, int a_n_cols)
    : VecFunctor(a_n_rows, a_n_cols), eqn(a_eqn), integ(a_integ) {}

    arma::mat operator()(double x, double y) const
    {
        return (eqn.*integ)(x, y);
    }
};

vec Equation::boundaryIntegral(int i, const MeshFn &u, int deg)
{
    int nv1, a1, b1;
    int j;
    vec integ = zeros<vec>(nc);
    
    nv1 = msh.p[i].size();
    // loop over all edges of the polygon we're in
    for (j = 0; j < nv1; j++)
    {
        a1 = msh.p[i][j];
        b1 = msh.p[i][(j+1)%nv1];
        
        
        iPlus = msh.p2p[i][j];
        if(iPlus >= 0)
        {
            UPlus = u.a.slice(iPlus);
        } else
        {
            // how to deal with BCs in a good way???
            //int p2 = dynamic_cast<PeriodicMesh &>(msh).bc[-msh.p2p[i][j]-1].p2;
            //boundaryTerm->UPlus = u.a.slice(p2);
        }
        
        msh.getOutwardNormal(i, a1, b1, nx, ny);
        integ += Quadrature::lineIntegral(msh, *boundaryTerm, a1, b1, deg*2);
    }
    
    return integ;
}

vec Equation::volumeIntegral(int i, int deg)
{
    return Quadrature::polygonIntegral(msh, *volumeTerm, i, deg*2);
}

MeshFn Equation::assemble(const MeshFn &u, double a_t/* = 0 */)
{
    int i, j;
    int deg = u.deg;
    int basisSize = (deg+1)*(deg+2)/2;
    double w, h;
    
    psi = zeros<vec>(basisSize);
    
    MeshFn b(msh, deg, nc);
    
    m = deg+1;
    
    // loop over all polygons
    for (i = 0; i < msh.np; i++)
    {
        iMinus = i;
        
        w = msh.bb[i][2];
        h = msh.bb[i][3];
        
        UMinus = u.a.slice(i);
        
        for (j = 0; j < basisSize; j++)
        {
            psi[j] = 1.0;
            
            psi_x = LegDerX(deg + 1, psi) * 2.0/w;
            psi_y = LegDerY(deg + 1, psi) * 2.0/h;
            
            b.a.slice(i).row(j) = volumeIntegral(i, deg).t();
            b.a.slice(i).row(j) -= boundaryIntegral(i, u, deg).t();
          
            psi[j] = 0.0;
        }
    }
    
    return b;
}

Jacobian Equation::jacobian(const MeshFn &f, double a_t)
{
    int i, j, k, e;
    int diagonalBlock, blockIdx;
    int deg = f.deg;
    int basisSize = (deg+1)*(deg+2)/2;
    int nv1, a1, b1;
    double w, h;
    uvec componentIndices = linspace<uvec>(0, (nc-1)*basisSize, nc);
    
    phi = zeros<vec>(basisSize);
    psi = zeros<vec>(basisSize);
    
    m = deg+1;
    
    Jacobian J(msh, f.deg, nc);
    
    for (i = 0; i < msh.np; i++)
    {
        w = msh.bb[i][2];
        h = msh.bb[i][3];
        diagonalBlock = J.rowBlock[i];
        
        U = f.a.slice(i);
        iPsi = i;
        
        for (j = 0; j < basisSize; j++)
        {
            phi(j) = 1;
            
            for (k = 0; k < basisSize; k++)
            {
                psi(k) = 1;
                psi_x = LegDerX(deg + 1, psi) * 2.0/w;
                psi_y = LegDerY(deg + 1, psi) * 2.0/h;
                
                J.blocks[diagonalBlock](componentIndices+k, componentIndices+j) 
                    += Quadrature::polygonIntegral(msh, *volumeJacobian, i, deg*2);
                
                nv1 = msh.p[i].size();
                
                blockIdx = diagonalBlock + 1;
                
                // loop over all edges of the polygon we're in
                for (e = 0; e < nv1; e++)
                {
                    a1 = msh.p[i][e];
                    b1 = msh.p[i][(e+1)%nv1];
                    
                    // get the outward facing normal vector
                    msh.getOutwardNormal(i, a1, b1, nx, ny);
                    
                    // get the neighbor corresponding to this edge
                    neighbor = msh.p2p[i][e];
                    iPhi = neighbor;
                    
                    // negative index indicates exterior edge
                    if (neighbor >= 0)
                    {
                        UNeighbor = f.a.slice(neighbor);
                        J.blocks[blockIdx](componentIndices+k, componentIndices+j) -= 
                            Quadrature::lineIntegral(msh, *boundaryJacobian, a1, b1, deg*2);
                        blockIdx++;
                    }
                    
                    iPhi = i;
                    J.blocks[diagonalBlock](componentIndices+k, componentIndices+j) -= 
                        Quadrature::lineIntegral(msh, *boundaryJacobian, a1, b1, deg*2);
                    
                }
                
                psi(k) = 0;
            }
            phi(j) = 0;
        }
    }
    
    return J;
}

Equation::Equation(PolyMesh &a_msh, int a_nc) : msh(a_msh), nc(a_nc)
{
    volumeTerm = new IntegrandFunctor(*this, &Equation::computeVolumeTerm, nc, 1);
    boundaryTerm = new IntegrandFunctor(*this, &Equation::computeBoundaryTerm, nc, 1);
    volumeJacobian = new IntegrandFunctor(*this, &Equation::computeVolumeJacobian, nc, nc);
    boundaryJacobian = new IntegrandFunctor(*this, &Equation::computeBoundaryJacobian, nc, nc);
}

Equation::~Equation()
{
    delete volumeTerm;
    delete boundaryTerm;
    delete volumeJacobian;
    delete boundaryJacobian;
}
