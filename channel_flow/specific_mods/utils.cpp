#include "tfa/Opers.h"
#include "tfa/Simulation.h"

#include <vector>
#include <fstream>

TF_Func bool stoppingCriterion(tfa::Simulation &sim)
{
  if(sim.IOParamD["_ElapsedTime"] > sim.IOParamD["_MaxTime"]) return tfa::Iter_Stop;
  
  return tfa::Iter_Continue;
}

TF_Func void init_props(tfa::Simulation &sim)
{
    double Re_tau = sim.IOParamD["Re_tau"];
    double &kinVisc = sim.IOParamD["kinVisc"];
    double &Ub = sim.IOParamD["Ub"];
    double &d = sim.IOParamD["delta"];
    double Ly = tfa::getLength(tfa::getSMesh(sim))[1];
    sim.IOParamD["_MaxTime"] = 1e5;
    sim.IOParamD["_RAn"] = Ub;
    sim.IOParamD["_RAn1"] = Ub;

    d = 0.5*Ly;
    Ub = d*0.5*Re_tau;
    kinVisc = 1.0/Re_tau;
    auto &msx = tfa::getField(sim, "momSrcx_C");
    double gp = sim.IOParamD["gradp"];
    tfa::oper_setConst(msx, gp);

    tfa::info("init_props completed.\n");
}

double my_rand0(double)
{
    return (1.0-2.0*drand48());
}

TF_Func void randomize_u(tfa::Simulation &sim)
{
    srand48(tfa::mpiRank());
    auto &uy = tfa::getField(sim, "uy_N");
    auto &uz = tfa::getField(sim, "uz_N");

    // Since we are passing 'my_rand0' as the argument to 'oper_apply', and
    // 'my_rand0' is a function that maps a double into a double (as opposed to
    // mapping a tfa::Vn into a tfa::Vn), it will work irrespective of how many
    // components the field has: it will be applied for every component of
    // every element.
    tfa::oper_apply(uy, my_rand0);
    tfa::oper_apply(uz, my_rand0);
    tfa::info("randomize_u completed\n");
}

TF_Func void init_profile_ux(tfa::Simulation &sim)
{
    auto Ub = sim.IOParamD["Ub"];
    auto d = sim.IOParamD["delta"];
    auto dim = tfa::getField(sim, "ux_N").dim;
    auto profile = [=](double, double y, double) -> tfa::Vn
    {
        // We want to apply the same profile for all components / simulations,
        // so we just replicate it 'dim' times.
        tfa::Vn result(dim, Ub*y/d*(2.0-y/d)+0.01*Ub*my_rand0(.1));
        return result;
    };
    tfa::initField(sim, "ux_N", profile);
    tfa::info("init_profile_ux completed.\n");

    auto &INF = tfa::getMatrix(sim,"Interp_NF");

    auto &ux = tfa::getField(sim,"ux_N");
    auto &uy = tfa::getField(sim,"uy_N");
    auto &uz = tfa::getField(sim,"uz_N");

    auto &uxf = tfa::getField(sim,"ux_F");
    auto &uyf = tfa::getField(sim,"uy_F");
    auto &uzf = tfa::getField(sim,"uz_F");

    tfa::oper_prod(INF,ux,uxf);
    tfa::oper_prod(INF,uy,uyf);
    tfa::oper_prod(INF,uz,uzf);
}

