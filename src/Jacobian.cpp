#include "Jacobian.h"

using namespace arma;

Jacobian::Jacobian(const PolyMesh &a_msh, int deg, int a_nc) : msh(a_msh)
{
    int i;
    int j;
    int numBlocks;
    
    n_rows = msh.np;
    nc = a_nc;
    basisSize = (deg + 1)*(deg + 2)/2;
    
    numBlocks = 0;
    
    for (i = 0; i < nb; i++)
    {
        blocks.push_back(mat());
        colIndices.push_back(i);
        rowBlock.push_back(numBlocks);
        numBlocks++;
        
        for (j = 0; j < msh.p2p[i].size(); j++)
        {
            blocks.push_back(mat());
            colIndices.push_back(msh.p2p[i][j]);
            numBlocks++;
        }
    }
    
    rowBlock.push_back(numBlocks);
}
