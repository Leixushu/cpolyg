#include "Advection.h"


MeshFn Advection::assemble(MeshFn &f)
{
    int i, j;
    int deg = f.deg;
    int basisSize = (deg+1)*(deg+2)/2;
    
    // only one component for advection equation
    MeshFn b(msh, deg, 1);
    
    // loop over all polygons
    for (i = 0; i < msh.np; i++)
    {
        for (j = 0; j < basisSize; j++)
        {
            
        }
    }
    
    return b;
    
//     b = np.zeros((msh.np * basis_size))
//     psi = np.zeros((basis_size))
//     psi_x = np.zeros((basis_size))
//     psi_y = np.zeros((basis_size))
//     
//     # loop over all polygons in the mesh
//     for i,p in enumerate(msh.p):
//         w = msh.bb[i,1,0]
//         h = msh.bb[i,1,1]
//         
//         # loop over all basis vectors
//         for j in range(basis_size):
//             # get the jth basis vector
//             psi[j] = 1
//             
//             # compute its gradient
//             fn.legderx(n_coeffs, psi, psi_x)
//             psi_x *= 2.0/w
//             
//             fn.legdery(n_coeffs, psi, psi_y)
//             psi_y *= 2.0/h
//             
//             fi = f[i,:]
//             bb = msh.bb[i]
//             v = msh.v[p]
//             
//             # set up product function
//             f_beta_dot_grad_psi = lambda x, y: (fn.eval_fn_p(fi, x, y, bb)
//                 *(beta(x,y)[0]*fn.eval_fn_p(psi_x, x, y, bb)
//                 + beta(x,y)[1]*fn.eval_fn_p(psi_y, x, y, bb)))
//             
//             f_beta_dot_grad_psi = fn.fn_callback(f_beta_dot_grad_psi)
//             
//             # compute the integrals on the interior and over the boundary (numerical flux)
//             b[i*basis_size + j] = fn.integrate_p(f_beta_dot_grad_psi, v, msh.tri[i])
//             b[i*basis_size + j] -= flux(f, psi, msh, i)
//             
//             # reset the basis vector to zero
//             psi[j] = 0
//     
//     return b
}
