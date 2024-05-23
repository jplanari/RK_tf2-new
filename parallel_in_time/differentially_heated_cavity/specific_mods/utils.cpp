#include "tfa/Opers.h"
#include "tfa/Simulation.h"

#include <vector>

TF_Func void init_props(tfa::Simulation &sim)
{
    double Ra = sim.IOParamD["Ra"]; //Rayleigh
    double Pr = sim.IOParamD["Pr"]; //Prandtl

    double &kinVisc = sim.IOParamD["kinVisc"];
    double &lambda = sim.IOParamD["lambda"];

    kinVisc = sqrt(Pr/Ra);
    lambda = 1/sqrt(Pr*Ra);
    
/*
    tfa::info("TAVG_start=%d\n",sim.IOParamI["_TAVG_Start"]);
    
    tfa::info("pre: max iters=%d\n",sim.cfgRT.maxIters);    
    sim.cfgRT.maxIters = sim.IOParamI["_TAVG_Start"] + (sim.cfgRT.maxIters - sim.IOParamI["_TAVG_Start"])/sim.IOParamI["nSims"];
    tfa::info("post: max iters=%d\n",sim.cfgRT.maxIters);    

    sim.IOParamD["_MaxTime"] = 100;

    if(sim.IOParamI["nSims"] > 1)
      tfa::info("Parallel-in-time, with %d rhs. Maximum number of iterations: %d\n",sim.IOParamI["nSims"],sim.cfgRT.maxIters);
*/
    tfa::info("init_props completed.\n");

    tfa::info("mpi_size is %d\n", tfa::mpiSize());
}

double my_rand0(double)
{
    return (1e-1-2e-1*drand48());
}

TF_Func void init_fields(tfa::Simulation &sim)
{
    srand48(tfa::mpiRank());
    auto &ux = tfa::getField(sim, "ux_N");
    auto &uy = tfa::getField(sim, "uy_N");
    auto &uz = tfa::getField(sim, "uz_N");

    // Since we are passing 'my_rand0' as the argument to 'oper_apply', and
    // 'my_rand0' is a function that maps a double into a double (as opposed to
    // mapping a tfa::Vn into a tfa::Vn), it will work irrespective of how many
    // components the field has: it will be applied for every component of
    // every element.
    tfa::oper_apply(ux, my_rand0);
    tfa::oper_apply(uy, my_rand0);
    tfa::oper_apply(uz, my_rand0);

    auto &INF = tfa::getMatrix(sim, "Interp_NF");
    auto &ufx = tfa::getField(sim, "ux_F");
    auto &ufy = tfa::getField(sim, "uy_F");
    auto &ufz = tfa::getField(sim, "uz_F");

    tfa::oper_prod(INF,ux,ufx);
    tfa::oper_prod(INF,uy,ufy);
    tfa::oper_prod(INF,uz,ufz);

}

TF_Func void init_omega(tfa::Simulation &sim)
{
    // For this particular application we need a field with the volume of each
    // cell, replicated in all the components. This volume field will be used
    // in the FSM, to scale the divergence of the predictor velocity. This is
    // needed because we use the laplacian multiplied by the volume as the
    // pressure equation, in order to make the matrix SPD.
    auto dim    = tfa::getField(sim, "ux_N").dim;
    auto &omega = tfa::getOrCreateField(sim, dim, "Omega_C", "Cells");
    const auto &m = tfa::getSMesh(sim);
    auto cells = tfa::getNumCells(m);
    std::vector<double> buffer(tfa::getNumEntries(omega));

    for (uint32_t c = 0; c < cells; c++)
    {
        const auto &C = tfa::getCellIJKFromId(m, c);
        double v = tfa::calcCellVolume(m, C);
        for (uint32_t d = 0; d < dim; d++)
        {
            buffer[c*dim+d] = v;
        }
    }
    tfa::oper_setData(omega, buffer.data());
    tfa::info("init_omega completed.\n");
}

TF_Func bool monitor(tfa::Simulation &sim)
{
    static bool first = true;
    static bool steady_state = false;
    const uint32_t ndim = tfa::getField(sim,"ux_N").dim;
    if (ndim>1){
    if (first)
    {
        tfa::info("# "
                  " ite wclock "
                  " ts time "
                  " resSolver "
                  " numPresIters"
                  "\n");
        first = false;
    }

    if (sim.IOParamI["_Iter"]%1 == 0)
    {
    static auto s = tfa::getSolver(sim, "Pressure_Solver");
    tfa::info("%d %.8e %.5e %.5e %.5e %d\n",
              sim.IOParamI["_Iter"], tfa::runTime(),
              sim.IOParamD["_TimeStep"], sim.IOParamD["_ElapsedTime"],
              tfa::oper_solve_residual(s), tfa::oper_solve_numIters(s));
    }
    }
    else
    {
    if (first)
    {
        tfa::info("# "
                  " ite wclock "
                  " ts time "
                  " resSolver "
		  " max(ux) "
      " numPresIters "
                  "\n");
        first = false;
    }
    double max_x;
    if (sim.IOParamI["_Iter"]%1 == 0)
    {
    auto &ux = tfa::getField(sim, "ux_N");
    max_x = std::max(fabs(tfa::oper_max(ux)[0]), fabs(tfa::oper_min(ux)[0]));
    static auto s = tfa::getSolver(sim, "Pressure_Solver");
    tfa::info("%d %.8e %.5e %.5e %.5e %.5e %d\n",
              sim.IOParamI["_Iter"], tfa::runTime(),
              sim.IOParamD["_TimeStep"], sim.IOParamD["_ElapsedTime"],
              tfa::oper_solve_residual(s),max_x, tfa::oper_solve_numIters(s));
    }
   
    }
    return tfa::Iter_Continue;
}

TF_Func void setupManualBocos(tfa::Simulation &sim)
{
auto &u = tfa::getField(sim, "ux_N");
auto cells = tfa::getDomainSize(sim, "Cells");
std::vector<double> buffer(tfa::getNumEntries(u), 1.0);
for (uint32_t it = cells*u.dim; it < buffer.size(); it++)
{
    buffer[it] = 0.0;
}
auto &cellMask = tfa::getOrCreateField(sim, "bnd", u);
tfa::oper_setData(cellMask, buffer.data());
}

TF_Func bool applyManualBocos(tfa::Simulation &sim)
{
    auto &bnd = tfa::getField(sim,"bnd");
	
    auto &ux = tfa::getField(sim,"ux_N");
    auto &uy = tfa::getField(sim,"uy_N");
    auto &uz = tfa::getField(sim,"uz_N");

    tfa::oper_prod(bnd,ux);
    tfa::oper_prod(bnd,uy);
    tfa::oper_prod(bnd,uz);

    return tfa::Iter_Continue;
}

TF_Func void calc_div_u(tfa::Simulation &sim)
{
    auto &DX = tfa::getMatrix(sim, "DivX_FC");
    auto &DY = tfa::getMatrix(sim, "DivY_FC");
    auto &DZ = tfa::getMatrix(sim, "DivZ_FC");

    auto &ux = tfa::getField(sim, "ux_F");
    auto &uy = tfa::getField(sim, "uy_F");
    auto &uz = tfa::getField(sim, "uz_F");

    auto &div = tfa::getOrCreateField(sim, "div(u)", DX, ux);
    tfa::oper_prod(DX, ux, div, 1.0, 0.0);
    tfa::oper_prod(DY, uy, div, 1.0, 1.0);
    tfa::oper_prod(DZ, uz, div, 1.0, 1.0);
}

TF_Func void computeTimes(tfa::Simulation &sim)
{
    CPUTimer_Info();
    CPUTimer_Reset();
}

TF_Func void startMainTimer(tfa::Simulation &sim)
{
    CPUTimer_Beg("main");
}

TF_Func void endMainTimer(tfa::Simulation &sim)
{
    CPUTimer_End("main");
}



