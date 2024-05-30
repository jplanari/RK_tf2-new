#include "tfa/Opers.h"
#include "tfa/Simulation.h"

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

TF_Func bool monitor(tfa::Simulation &sim)
{
    static bool first = true;
    static bool steady_state = false;
    const uint32_t ndim = tfa::getField(sim,"ux").dim;
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
    auto &ux = tfa::getField(sim, "ux");
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

TF_Func void calc_div_u(tfa::Simulation &sim)
{
    auto &DX = tfa::getMatrix(sim, "DivX_FN");
    auto &DY = tfa::getMatrix(sim, "DivY_FN");
    auto &DZ = tfa::getMatrix(sim, "DivZ_FN");

    auto &ux = tfa::getField(sim, "ux_F");
    auto &uy = tfa::getField(sim, "uy_F");
    auto &uz = tfa::getField(sim, "uz_F");

    auto &div = tfa::getOrCreateField(sim, "div(u)", DX, ux);
    tfa::oper_prod(DX, ux, div, 1.0, 0.0);
    tfa::oper_prod(DY, uy, div, 1.0, 1.0);
    tfa::oper_prod(DZ, uz, div, 1.0, 1.0);
}

void printPhi(tfa::Simulation &sim)
{
  std::string Re = std::to_string(sim.IOParamD["Re_tau"]);
  std::string Rkfct = std::to_string(sim.IOParamD["RKfct"]);
  std::string RKmethod = sim.IOParamS["RKmethod"];

  std::string filename = "results/phi_Re_tau-"+Re+"_RKfct-"+Rkfct+"_"+RKmethod+".dat";
  std::ofstream file;

  static bool first = true;
  if (first) {
    file.open(filename);
    first = false;
  }
  else
    file.open(filename,std::ios::app);
  
  if (file.is_open()){
    file << sim.IOParamI["_Iter"] << "\t" << atan(sim.IOParamD["_EVimag"]/sim.IOParamD["_EVreal"]) << std::endl; 
    file.close();;
  }
}

