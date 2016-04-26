#include "TimeIntegration.h"
#include "Preconditioners.h"

#define kNewtonMaxIterations 100
#define kNewtonTolerance 5.e-13

using namespace std;
using namespace arma;

MeshFn ForwardEuler::advance(const MeshFn &u, const double dt, const double t) {
    return u + dt * M.solve(eqn.assemble(u, t));
}

MeshFn RK4::advance(const MeshFn &u, const double dt, const double t) {
    CH_TIMERS("RK4Advance");
    MeshFn k1 = dt*M.solve(eqn.assemble(u, t));
    MeshFn k2 = dt*M.solve(eqn.assemble(u + 0.5*k1, t + 0.5*dt));
    MeshFn k3 = dt*M.solve(eqn.assemble(u + 0.5*k2, t + 0.5*dt));
    MeshFn k4 = dt*M.solve(eqn.assemble(u + k3, t + dt));
    
    return u + (k1 + 2*k2 + 2*k3 + k4)/6.0;
}

MeshFn RK2::advance(const MeshFn &u, const double dt, const double t) {
    MeshFn k1 = dt*M.solve(eqn.assemble(u, t));
    MeshFn k2 = dt*M.solve(eqn.assemble(u + 0.5*k1, t + 0.5*dt));
    
    return u + 0.5*(k1 + k2);
}

MeshFn BackwardEuler::advance(const MeshFn &u, const double dt,
                              const double t) {
    MeshFn unp1 = u;
    MeshFn r;
    MeshFn b;
    int k;

    // Newton solve
    cout << "  Beginning Newton solve" << endl;
    for (k = 0; k < kNewtonMaxIterations; k++) {

        b = eqn.assemble(unp1, t + dt);
        r = M.matvec(unp1 - u) - dt * b;

        b.gnuplot("plt/b" + to_string(k) + ".gnu");

        cout << "    Iteration number " << k
             << ", residual norm = " << r.L2Norm().max() << endl;

        if (r.L2Norm().max() < kNewtonTolerance) {
            cout << "    Converged to tolerance" << endl;
            break;
        }

        Jacobian B = eqn.jacobian(unp1, t + dt);
        B *= -dt;
        B += M;
        BlockILU0 pc(B);
        //BlockJacobi pc(B);

        unp1 -= B.solve(r, pc, kGMRESSolver);
        //unp1 -= B.solve(r, pc, kJacobiSolver);
    }

    return unp1;
}

void getDIRK(int nStages, arma::mat &A, arma::vec &b, arma::vec &c) {
    static double const alpha=0.435866521508459;
    switch(nStages) {
        case 3:
            A = {{(1+alpha)/2.0, 0,             0},
                {-alpha/2.0,    (1+alpha)/2.0, 0},
                {1+alpha,      -(1+2*alpha),   (1+alpha)/2.0}};
            b = {1.0/6.0/(alpha*alpha),1-1.0/3.0/(alpha*alpha), 1.0/6.0/(alpha*alpha)};
            c = {(1+alpha)/2.0, 1.0/2.0, (1-alpha)/2.0};
            break;
        case 5:
            // Cooper and Sayfy (1979)
            double g = (6-sqrt(6))/10.0;
            A = {{g,0,0,0,0},
                 {(-6+5*sqrt(6))/14.0,g,0,0,0},
                 {(888+607*sqrt(6))/2850.0,(126-161*sqrt(6))/1425.0,g,0,0},
                 {(3153-3082*sqrt(6))/14250.0,(3213+1148*sqrt(6))/28500.0,(-267+88*sqrt(6))/500.0,g,0},
                 {(-32583+14638*sqrt(6))/71250.0,(-17199+364*sqrt(6))/142500.0,(1329-544*sqrt(6))/2500.0,(-96+131*sqrt(6))/625.0,g}};
            b = {0, 0, 1/9.0, (16 - sqrt(6))/36.0, (16 + sqrt(6))/36.0};
            c = {g, (6+9*sqrt(6))/35.0, 1, (4-sqrt(6))/10.0, (4+sqrt(6))/10.0};
            
            break;
    }
}

MeshFn DIRK::advance(const MeshFn &u, const double dt, const double t)
{
    CH_TIMERS("DIRK");
    MeshFn zero(*u.msh, u.deg, u.nc);
    zero.a.fill(0);
    //array<MeshFn, nStages> k = {{zero, zero, zero}};
    field<MeshFn> k(nStages);
    k.fill(zero);
    
    for (int stage = 0; stage < nStages; stage++)
    {
        cout << "  DIRK stage " << stage+1 << endl;
        
        cout << "  Beginning Newton solve" << endl;
        for (int iter = 0; iter < kNewtonMaxIterations; iter++)
        {
            MeshFn cu = u;
            for (int j = 0; j <= stage; j++)
            {
                cu += A(stage, j)*k(j);
            }
            
            MeshFn rhs = eqn.assemble(cu, t + dt*c(stage));
            MeshFn r = M.matvec(k(stage)) - dt*rhs;
            
            cout << "    Iteration number " << iter << ", residual norm = "
                 << r.L2Norm().max() << endl;
            
            if (r.L2Norm().max() < kNewtonTolerance)
            {
                cout << "    Converged to tolerance" << endl;
                break;
            }
            
            Jacobian B = eqn.jacobian(cu, t + dt*c(stage));
            B *= -dt*A(stage, stage);
            B += M;
            
            CH_TIMER("Linear solve", t1);
            CH_START(t1);
            
            BlockILU0 pc(B);
            
            //NoPreconditioner pc;
            k(stage) -= B.solve(r, pc, kGMRESSolver);
            
            CH_STOP(t1);
        }
    }
    
    MeshFn unp1 = u;
    
    for (int i = 0; i < nStages; i++)
    {
        unp1 += b(i)*k(i);
    }
    
    return unp1;
}

void getIRK(int nStages, arma::mat &A, arma::mat &Ainv, arma::vec &b, arma::vec &c) {
    switch(nStages) {
        case 2:
            A = {{5/12.0, -1/12.0},
                 {3/4.0,   1/4.0}};
            Ainv = {{ 3/2.0, 1/2.0},
                    {-9/2.0, 5/2.0}};
            b = {3/4.0, 1/4.0};
            c = {1/3.0, 1};
            break;
        case 3:
            A = {{(88-7*sqrt(6))/360.0,
                    (296-169*sqrt(6))/1800.0,
                    (-2+3*sqrt(6))/225.0},
               {(296+169*sqrt(6))/1800.0,
                    (88+7*sqrt(6))/360.0,
                    (-2-3*sqrt(6))/225.0},
               {(16-sqrt(6))/36.0,
                    (16+sqrt(6))/36.0,
                    1/9.0}};
            Ainv = {{60*(1/30.0 + 1.0/(20*sqrt(6))),
                        60*(-1/50.0 + 29.0/(300*sqrt(6))),
                        60*(1/150.0 - sqrt(2/3.0)/75.0)},
                  {60*(-1/50.0 - 29.0/(300*sqrt(6))),
                        60*(1/30.0 - 1.0/(20*sqrt(6))),
                        60.0*(1/150.0 + sqrt(2/3.0)/75.0)},
                  {60*(-1/60.0 + (2*sqrt(2/3.0))/15.0),
                        60*(-1/60.0 - (2*sqrt(2/3.0))/15.0),
                        5}};
            b = {(16-sqrt(6))/36.0, (16+sqrt(6))/36.0, 1/9.0};
            c = {(4-sqrt(6))/10.0, (4+sqrt(6))/10.0, 1};
            break;
    }
}

MeshFn IRK::advance(const MeshFn &u, const double dt, const double t)
{
    CH_TIMERS("Full IRK - Regular Butcher Tableau");
    //const static int nStages = 3;
    MeshFn zero(*u.msh, u.deg, u.nc);
    zero.a.fill(0);
    
    field<MeshFn> k(nStages);
    field<MeshFn> r(nStages);
    field<MeshFn> w(nStages);
    
    k.fill(zero);
    
    field<Jacobian> Js(nStages);
    
    int vecSize = u.nc*(u.deg+1)*(u.deg+2)*0.5*u.msh->np;
    
    vec bigK = zeros(vecSize*nStages);
    vec bigR = zeros(vecSize*nStages);
    
    cout << "  Fully implicit IRK solve" << endl;
    cout << "  Beginning Newton solve" << endl;
    for (int iter = 0; iter < kNewtonMaxIterations; iter++) {
        double rNorm = 0;
        for (int stage = 0; stage < nStages; stage++) {
            w(stage) = u;
            for (int j = 0; j < nStages; j++) {
                w(stage) += A(stage, j)*k(j);
            }
            MeshFn rhs = eqn.assemble(w(stage), t + dt*c(stage));
            r(stage) = M.matvec(k(stage)) - dt*rhs;
            rNorm += pow(r(stage).L2Norm().max(), 2);
        }
        rNorm = sqrt(rNorm);
        cout << "    Iteration number " << iter << ", residual norm = "
             << rNorm << endl;
        if (rNorm < kNewtonTolerance) {
            cout << "    Converged to tolerance" << endl;
            break;
        }
        
        for (int stage = 0; stage < nStages; stage++) {
            Js(stage) = eqn.jacobian(w(stage), t + dt*c(stage));
        }
        
        field<BlockMatrix> JBlocks(nStages, nStages);
        for (int i = 0; i < nStages; i++) {
            Js(i) *= -dt;
            for (int j = 0; j < nStages; j++) {
                if (i != j) {
                    JBlocks(i,j) = Js(i);
                    JBlocks(i,j) *= A(i,j);
                }
            }
            Js(i) *= A(i,i);
            Js(i) += M;
            JBlocks(i,i) = Js(i);
        }
        BlockMatrix fullJ = BlockMatrix::blockBlockMatrix(JBlocks);
        
        for (int stage = 0; stage < nStages; stage++) {
            bigR.rows(stage*vecSize, (stage+1)*vecSize-1) = vectorise(r(stage).a);
        }
        
        CH_TIMER("Linear solve", t1);
        CH_START(t1);
        
        BlockILU0 pc(fullJ);
        //NoPreconditioner pc;
        
        double tol = 1.e-14;
        int maxIt = 20000;
        fullJ.gmres(bigR, bigK, 20, tol, maxIt, pc);
        
        CH_STOP(t1);
        
        cout << "    GMRES iterations: " << maxIt << endl;
        cube xCube(vecSize, 1, 1);
        
        for (int stage = 0; stage < nStages; stage++) {
            xCube.slice(0) = bigK.rows(stage*vecSize, (stage+1)*vecSize-1);
            
            k(stage).a -= reshape(xCube, (u.deg+1)*(u.deg+2)*0.5, u.nc, u.msh->np);
        }
    }
    
    MeshFn unp1 = u;
    for (int i = 0; i < nStages; i++)
    {
        unp1 += b(i)*k(i);
    }
    
    return unp1;
}

MeshFn IRK::newAdvance(const MeshFn &u, const double dt, const double t)
{
    CH_TIMERS("New IRK - Inverted Butcher Tableau");
    MeshFn zero(*u.msh, u.deg, u.nc);
    zero.a.fill(0);
    
    field<MeshFn> k = {zero, zero, zero};
    field<MeshFn> r = {zero, zero, zero};
    field<MeshFn> w = {zero, zero, zero};
    
    field<Jacobian> Js(nStages);
    
    int vecSize = u.nc*(u.deg+1)*(u.deg+2)*0.5*u.msh->np;
    
    vec bigU = zeros(vecSize*nStages);
    vec bigR = zeros(vecSize*nStages);
    
    cout << "  Fully implicit (inverted) IRK solve" << endl;
    cout << "  Beginning Newton solve" << endl;
    for (int iter = 0; iter < kNewtonMaxIterations; iter++) {
        for (int i = 0; i < nStages; i++) {
            k(i).a.fill(0);
            for (int j = 0; j < nStages; j++) {
                k(i) += Ainv(i, j)*w(j);
            }
        }
        
        double rNorm = 0;
        for (int stage = 0; stage < nStages; stage++) {
            MeshFn rhs = eqn.assemble(u + w(stage), t + dt*c(stage));
            r(stage) = M.matvec(k(stage)) - dt*rhs;
            rNorm += pow(r(stage).L2Norm().max(), 2);
        }
        rNorm = sqrt(rNorm);
        cout << "    Iteration number " << iter << ", residual norm = "
             << rNorm << endl;
        if (rNorm < kNewtonTolerance) {
            cout << "    Converged to tolerance" << endl;
            break;
        }
        
        for (int stage = 0; stage < nStages; stage++) {
            Js(stage) = eqn.jacobian(u + w(stage), t + dt*c(stage));
        }
        
        field<BlockMatrix> JBlocks(nStages, nStages);
        for (int i = 0; i < nStages; i++) {
            for (int j = 0; j < nStages; j++) {
                if (i != j) {
                    JBlocks(i,j) = M;
                    JBlocks(i,j) *= Ainv(i,j);
                }
            }
            
            MassMatrix Mscaled = M;
            Mscaled *= Ainv(i,i);
            Js(i) *= -dt;
            Js(i) += Mscaled;
            JBlocks(i,i) = Js(i);
        }
        BlockMatrix fullJ = BlockMatrix::blockBlockMatrix(JBlocks);
        
        field<BlockMatrix> JPCBlocks(nStages, nStages);
        for (int i = 0; i < nStages; i++) {
            for (int j = 0; j < nStages; j++) {
                if (i != j) {
                    JPCBlocks(i,j) = BlockMatrix(Js(i).bl, Js(i).n_rows);
                    //JPCBlocks(i,j) = M;
                    //JPCBlocks(i,j) *= Ainv(i, j);
                }
            }
            JPCBlocks(i,i) = Js(i);
        }
        BlockMatrix fullJPC = BlockMatrix::blockBlockMatrix(JPCBlocks);
        
        
        for (int stage = 0; stage < nStages; stage++) {
            bigR.rows(stage*vecSize, (stage+1)*vecSize-1) = vectorise(r(stage).a);
        }
        
        CH_TIMER("Linear solve", t1);
        CH_START(t1);
        
        //NoPreconditioner pc;
        BlockILU0 pc(fullJPC);
        //BlockILU0 pc(fullJ);
        
        //pc.DD.spy("plt/DD.gnu");
        //pc.DU.spy("plt/DU.gnu");
        //pc.L.spy("plt/L.gnu");
        
        double tol = 1.e-14;
        int maxIt = 20000;
        fullJ.gmres(bigR, bigU, 20, tol, maxIt, pc);
        
        CH_STOP(t1);
        
        cout << "    GMRES iterations: " << maxIt << endl;
        cube xCube(vecSize, 1, 1);
        
        for (int stage = 0; stage < nStages; stage++) {
            xCube.slice(0) = bigU.rows(stage*vecSize, (stage+1)*vecSize-1);
            
            w(stage).a -= reshape(xCube, (u.deg+1)*(u.deg+2)*0.5, u.nc, u.msh->np);
        }
    }
    
    mat btAinv = b.t()*Ainv;
    MeshFn unp1 = u;
    for (int i = 0; i < nStages; i++) {
        unp1 += btAinv(i)*(w(i));
    }
    
    return unp1;
}
