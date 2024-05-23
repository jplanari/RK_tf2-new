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

bool determineSteadyState(double max_x, double max_u, tfa::Simulation &sim)
{
  if(sim.IOParamI["_Iter"]%1000 == 0)
  {
      if(fabs(max_x-max_u)/max_u < 1e-1){
        sim.IOParamD["_MaxTime"] = sim.IOParamD["_ElapsedTime"] + sim.IOParamD["nFT"]*4*M_PI/(0.65*sim.IOParamD["maxU"]);
        sim.IOParamI["_TAVG_Start"] = sim.IOParamI["_Iter"];
        return true;
      }
      sim.IOParamD["maxU"] = max_x;
      return false;
  } 
}

bool rollingSteadyState(double max_x, int window_size, tfa::Simulation &sim)
{
  double curr_ra = sim.IOParamD["_RAn"];
  curr_ra += max_x;
  if(sim.IOParamI["_Iter"]%window_size==0)
  {
    curr_ra /= (double) window_size;
    if (fabs(curr_ra-sim.IOParamD["_RAn1"])/sim.IOParamD["_RAn1"] < 1e-1){
        sim.IOParamD["_MaxTime"] = sim.IOParamD["_ElapsedTime"] + sim.IOParamD["nFT"]*4*M_PI/(0.65*sim.IOParamD["maxU"]);
        sim.IOParamI["_TAVG_Start"] = sim.IOParamI["_Iter"];
        return true;
    }
    else
      sim.IOParamD["_RAn1"] = sim.IOParamD["_RAn"];
  }
  return false;
}

TF_Func bool monitor(tfa::Simulation &sim)
{
    static bool first = true;
    const uint32_t ndim = tfa::getField(sim,"ux_N").dim;
    if (ndim>1){
    if (first)
    {
        tfa::info("# "
                  " ite wclock "
                  " ts time "
                  " resSolver "
                  "\n");
        first = false;
    }

    if (sim.IOParamI["_Iter"]%10 == 0)
    {
    static auto s = tfa::getSolver(sim, "Pressure_Solver");
    tfa::info("real=%e, imag=%e\n",sim.IOParamD["_EVreal"],sim.IOParamD["_EVimag"]);
    tfa::info("%d %.8e %.5e %.5e %.5e\n",
              sim.IOParamI["_Iter"], tfa::runTime(),
              sim.IOParamD["_TimeStep"], sim.IOParamD["_ElapsedTime"],
              tfa::oper_solve_residual(s));
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
                  "\n");
        first = false;
    }
    double max_x;
    if (sim.IOParamI["_Iter"]%10 == 0)
    {
    printPhi(sim);
    auto &ux = tfa::getField(sim, "ux_N");
    max_x = std::max(fabs(tfa::oper_max(ux)[0]), fabs(tfa::oper_min(ux)[0]));
    static auto s = tfa::getSolver(sim, "Pressure_Solver");
    tfa::info("%d %.8e %.5e %.5e %.5e %.5e\n",
              sim.IOParamI["_Iter"], tfa::runTime(),
              sim.IOParamD["_TimeStep"], sim.IOParamD["_ElapsedTime"],
              tfa::oper_solve_residual(s),max_x);
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


