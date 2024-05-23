#include "tfa/Opers.h"
#include "tfa/Simulation.h"

static void init_cfl(tfa::Simulation &sim)
{
    double DX = tfa::meshCharLength(sim);
    double nu  = tfa::getParamD(sim, "kinVisc");
    double Pr  = tfa::getParamD(sim, "Pr");
    double cfl = tfa::getParamD(sim, "cfl");
    tfa::setParamD(sim, "dx", DX);
    tfa::setParamD(sim, "tvisc", 0.2*std::min(Pr, 1.0)*DX*DX/nu);
    tfa::setDt(sim, cfl*tfa::getParamD(sim, "tvisc"));
    tfa::info("timestep = %.5e\n", tfa::getParamD(sim, "_TimeStep"));
    tfa::info("tvisc = %.5e\n", tfa::getParamD(sim, "tvisc"));
}

TF_Func void init_const(tfa::Simulation &sim)
{
    double Ra = tfa::getParamD(sim, "Ra");
    double Pr = tfa::getParamD(sim, "Pr");
    double Tmin = tfa::getParamD(sim, "Tmin");
    double Tmax = tfa::getParamD(sim, "Tmax");
    double gx = tfa::getParamD(sim, "gx");
    double gy = tfa::getParamD(sim, "gy");
    double gz = tfa::getParamD(sim, "gz");

    double g = sqrt(gx*gx+gy*gy+gz*gz);
    double beta = 1.0/Tmin;
    double DT = Tmax-Tmin;
    double L = 1.0;
    double rho = 1.0;
    double k = 1.0;
    double kinVisc = sqrt(Pr*g*beta*DT*L*L*L/Ra);
    double cp = k*Pr/kinVisc/rho;

    tfa::info("Ra = %.5e\n", Pr*L*L*L*g*beta*DT/kinVisc/kinVisc);
    tfa::info("Pr = %.5e\n", rho*kinVisc*cp/k);
    tfa::info("cp = %.5e\n", cp);
    tfa::info("nu = %.5e\n", kinVisc);

    tfa::setParamD(sim, "beta",      beta);
    tfa::setParamD(sim, "kinVisc",   kinVisc);
    tfa::setParamD(sim, "thermCond", k);
    tfa::setParamD(sim, "rhoCp",     rho*cp);
    tfa::setParamD(sim, "Tref",      Tmin);

    init_cfl(sim);
}

TF_Func bool cfl_condition(tfa::Simulation &sim)
{
    auto &ux = tfa::getField(sim, "ux");
    auto &uy = tfa::getField(sim, "uy");
    auto &uz = tfa::getField(sim, "uz");
    auto modv = tfa::getTmpFieldS(sim, ux);

    tfa::oper_fmadd(ux, ux, *modv, 1.0, 0.0);
    tfa::oper_fmadd(uy, uy, *modv, 1.0, 1.0);
    tfa::oper_fmadd(uz, uz, *modv, 1.0, 1.0);

    double uref  = sqrt(tfa::max(tfa::oper_max(*modv)));
    double tvisc = tfa::getParamD(sim, "tvisc");
    double DX    = tfa::getParamD(sim, "dx");
    double cfl   = tfa::getParamD(sim, "cfl");
    tfa::setDt(sim, cfl*std::min(tvisc, 0.35*DX/uref));

    return tfa::Iter_Continue;
}

TF_Func bool energy_const(tfa::Simulation &sim)
{
    auto &T   = tfa::getField(sim, "T");
    auto &R   = tfa::getOrCreateField(sim, "R", T);
    auto &R0  = tfa::getOrCreateField(sim, "R0", T);
    auto &T0  = tfa::getOrCreateField(sim, "T0", T);
    auto &ufx = tfa::getField(sim, "ux_F");
    auto &ufy = tfa::getField(sim, "uy_F");
    auto &ufz = tfa::getField(sim, "uz_F");

    // -------------------------------------------------------------------------

    auto &L = tfa::getMatrix(sim, "DiffT_NN");

    double thermCond = tfa::getParamD(sim, "thermCond");
    double rhoCp     = tfa::getParamD(sim, "rhoCp");

    tfa::oper_prod(L, T, R, thermCond/rhoCp);

    // -------------------------------------------------------------------------

    auto &DX  = tfa::getMatrix(sim, "DivX_FN");
    auto &DY  = tfa::getMatrix(sim, "DivY_FN");
    auto &DZ  = tfa::getMatrix(sim, "DivZ_FN");
    auto &INF = tfa::getMatrix(sim, "InterpSP_NF");

    auto TF = tfa::getTmpFieldS(sim, ufx);
    auto xx = tfa::getTmpFieldS(sim, ufx);

    tfa::oper_prod(INF, T, *TF);

    tfa::oper_fmadd(*TF, ufx, *xx, 1.0, 0.0);
    tfa::oper_prod(DX, *xx, R, -1.0, 1.0);
    tfa::oper_fmadd(*TF, ufy, *xx, 1.0, 0.0);
    tfa::oper_prod(DY, *xx, R, -1.0, 1.0);
    tfa::oper_fmadd(*TF, ufz, *xx, 1.0, 0.0);
    tfa::oper_prod(DZ, *xx, R, -1.0, 1.0);

    // -------------------------------------------------------------------------

    double ts = tfa::getDt(sim);
    tfa::oper_axpy(R,  T, +1.5*ts, 1.0);
    tfa::oper_axpy(R0, T, -0.5*ts, 1.0);

    tfa::oper_axpy(T, T0, 1.0, -1.0);
    tfa::setParamD(sim, "TEQ_Delta_T", sqrt(tfa::max(tfa::oper_dot(T0, T0))));

    tfa::oper_copy(R, R0);
    tfa::oper_copy(T, T0);

    return tfa::Iter_Continue;
}

TF_Func bool couple_energy_momentum(tfa::Simulation &sim)
{
    double gx   = tfa::getParamD(sim, "gx");
    double gy   = tfa::getParamD(sim, "gy");
    double gz   = tfa::getParamD(sim, "gz");
    double beta = tfa::getParamD(sim, "beta");
    double Tref = tfa::getParamD(sim, "Tref");

    auto &srcx = tfa::getField(sim, "momSrcx");
    auto &srcy = tfa::getField(sim, "momSrcy");
    auto &srcz = tfa::getField(sim, "momSrcz");
    auto &T    = tfa::getField(sim, "T");

    auto Tx = tfa::getTmpFieldS(sim, T);
    tfa::oper_copy(T, *Tx);

    // Tx = Tref-TC
    tfa::oper_add(*Tx, -Tref);

    // S[i] = g[i]*beta*(Tref-TC)
    tfa::oper_axpy(*Tx, srcx, -gx*beta, 0.0);
    tfa::oper_axpy(*Tx, srcy, -gy*beta, 0.0);
    tfa::oper_axpy(*Tx, srcz, -gz*beta, 0.0);

    return tfa::Iter_Continue;
}

TF_Func bool decider_nf(tfa::Simulation &sim)
{
    double runTime = tfa::runTime();
    static auto S = tfa::getSolver(sim, "Pressure_Solver");
    tfa::info("%d %.8f resP:%.5e (%3d) DeltaT:%.5e ts:%.5e et:%.8e\n"
              , tfa::getIter(sim), runTime
              , tfa::oper_solve_residual(S)
              , tfa::oper_solve_numIters(S)
              , tfa::getParamD(sim, "TEQ_Delta_T")
              , tfa::getDt(sim)
              , tfa::getElapsedTime(sim)
              );
    return tfa::Iter_Continue;
}
