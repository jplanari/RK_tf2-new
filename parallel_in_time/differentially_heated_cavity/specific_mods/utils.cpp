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

