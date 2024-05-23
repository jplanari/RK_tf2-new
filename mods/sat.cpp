#include "tfa/Simulation.h"
#include "tfa/Opers.h"
#include <cmath>

#include "fsm_rk/butcherTableaux.h"
#include "sat/stability.h"
#include "sat/eigenbounds.h"
#include "fsm_rk/multiRK.h"

void computeDT(tfa::Simulation &sim)
{
    double ev_r = sim.IOParamD["_EVreal"];
    double ev_i = sim.IOParamD["_EVimag"];
    double ev_norm = sqrt(ev_r*ev_r+ev_i*ev_i);

    double phi = M_PI-atan(ev_i/ev_r);

    sim.IOParamD["_TimeStep"] = sim.IOParamD["RKfct"]*stabilityRegion(phi,intScheme(sim.IOParamS["RKmethod"]))/ev_norm;
}

void computeDT_efficiency(tfa::Simulation &sim)
{
    double ev_r = sim.IOParamD["_EVreal"];
    double ev_i = sim.IOParamD["_EVimag"];
    double ev_norm = sqrt(ev_r*ev_r+ev_i*ev_i);

    if(sim.IOParamS["RKmethod"] != "efficiency") computeDT(sim);
    else{
      long unsigned int idd=0;
      double phi = M_PI-atan(ev_i/ev_r);
      double h;
      double hmax=0.0, htilde, hh=0.0;
      butcherTableau candidate;
      double f;
    
      if(sim.IOParamS["FPJ"] == "yes") f = 0.85;
      else f = 0.0;

      std::vector<std::string> pool = multiRK_method(sim);
	    
      for(long unsigned int id=0; id<pool.size(); ++id){
        candidate = intScheme(pool.at(id));
	      h = stabilityRegion(phi,candidate);
	      htilde = (h/ev_norm)/tauMethod(candidate.b.size(),f);
	      if(htilde>hmax){
		      hmax = htilde;
          hh = h/ev_norm;
		      idd = id;
	      }
	    }
      
      tfa::info("Scheme selected = %d\n", idd);
      sim.IOParamI["_schEff"] = (int)idd; //idea is to remove this option, kept now just to check if everything works fine
      sim.IOParamS["RKname"] = pool.at(idd);
      sim.IOParamD["_TimeStep"] = sim.IOParamD["RKfct"]*hh;
    }

    if(sim.IOParamS["energy"] == "yes"){
      double ev_i_T = ev_i/sim.IOParamD["Pr"];
      double ev_norm_T = sqrt(ev_r*ev_r+ev_i_T*ev_i_T);
      butcherTableau energyTab;
      if(sim.IOParamS["RKmethod"] == "efficiency"){
        std::vector<std::string> pool = multiRK_method(sim);
        energyTab = intScheme(pool.at((long unsigned int)sim.IOParamI["_schEff"]));
      }
      else energyTab = intScheme(sim.IOParamS["RKmethod"]);
  
      double phi_T = M_PI-atan(ev_i_T/ev_r);
  
      sim.IOParamD["_TimeStep"] = std::min(sim.IOParamD["_TimeStep"],sim.IOParamD["RKfct"]*stabilityRegion(phi_T,intScheme(sim.IOParamS["RKmethod"]))/ev_norm_T);
    }
}

TF_Func void SetUp_SAT_Gershgorin_efficiency(tfa::Simulation &sim)
{
  computeEV_Gershgorin(sim);
	computeDT_efficiency(sim);
}

TF_Func bool Iter_SAT_Gershgorin_efficiency(tfa::Simulation &sim)
{
	computeImagEV_Gershgorin(sim);
	computeDT_efficiency(sim);
	return tfa::Iter_Continue;
}

TF_Func void SetUp_SAT_GershgorinMat_efficiency(tfa::Simulation &sim)
{
  computeRealEV_GershgorinMat(sim);
  computeImagEV_GershgorinMat(sim);
	computeDT_efficiency(sim);
}

TF_Func bool Iter_SAT_GershgorinMat_efficiency(tfa::Simulation &sim)
{
	computeImagEV_GershgorinMat(sim);
	computeDT_efficiency(sim);
	return tfa::Iter_Continue;
}
TF_Func void SetUp_SAT_AlgEigCD_efficiency(tfa::Simulation &sim)
{
  computeEV_AlgEigCD(sim);
	computeDT_efficiency(sim);
}

TF_Func bool Iter_SAT_AlgEigCD_efficiency(tfa::Simulation &sim)
{
	computeImagEV_AlgEigCD(sim);
	computeDT_efficiency(sim);
	return tfa::Iter_Continue;
}
