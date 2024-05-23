#include "tfa/ll_Opers.h"
#include "tfa/Opers.h"
#include "tfa/Simulation.h"

#include "fsm_rk/butcherTableaux.h"
#include "fsm_rk/setup_fsm.h"
#include "fsm_rk/cd.h"
#include "fsm_rk/fsm_rk.h"

#include "energy/energy.h"

#include <vector>
#include <cstdarg>

// +---------------------------------------------------------------------------+
// |                              General remarks                              |
// +---------------------------------------------------------------------------+

// This module introduces some functions to help solving the Navier-Stokes
// equations. For the moment only incompressible fluids have been considered.

TF_Func void SetUp_Momentum_RK(tfa::Simulation &sim)
{
    // This function exists basically to initialize all the fields required by
    // the official FSM, with the correct number of components, i.e., one
    // component per simulation.
    TF_uAssert((sim.meshes.size() == 1) xor tfa::hasSMesh(sim), "not yet implemented for multiple meshes");
    TF_uAssert(tfa::hasField(sim, "ux_N"), "required field ux_N not defined");
    TF_uAssert(tfa::hasField(sim, "uy_N"), "required field uy_N not defined");
    TF_uAssert(tfa::hasField(sim, "uz_N"), "required field uz_N not defined");
    TF_uAssert(tfa::hasField(sim, "Omega_C"), "required field Omega_C not defined");

    const auto &meshName = tfa::getMeshName(sim);
    const auto &msuf     = tfa::meshSuffix(meshName);
    const auto &cells    = "Cells"+msuf;
    const auto &faces    = "Faces"+msuf;
    const auto &nodes    = "Nodes"+msuf;

    const auto smesh = tfa::hasSMesh(sim);
   
    const auto &pMod = smesh? "smeshPatterns.so" : "patterns.so";
    const auto &kMod = smesh? "smeshInterpolators.so" : "interpolators.so";
    
    if (not tfa::hasMatrix(sim, "Id_NC"))
    {
        tfa::newMatrix(sim, "Id_NC", "pId_NC", pMod, "kId_NC", kMod, meshName);
    }
    if (not tfa::hasMatrix(sim, "Id_CN"))
    {
        tfa::newMatrix(sim, "Id_CN", "pId_CN", pMod, "kId_CN", kMod, meshName);
    }

    const auto dim = tfa::getField(sim, "ux_N").dim;

    tfa::getOrCreateField(sim, dim, "upredx_C",  cells);
    tfa::getOrCreateField(sim, dim, "upredy_C",  cells);
    tfa::getOrCreateField(sim, dim, "upredz_C",  cells);

    tfa::getOrCreateField(sim, dim, "momSrcx_C", cells);
    tfa::getOrCreateField(sim, dim, "momSrcy_C", cells);
    tfa::getOrCreateField(sim, dim, "momSrcz_C", cells);

    tfa::getOrCreateField(sim, dim, "P_C",       cells);

    tfa::getOrCreateField(sim, dim, "ux_F", faces);
    tfa::getOrCreateField(sim, dim, "uy_F", faces);
    tfa::getOrCreateField(sim, dim, "uz_F", faces);

    tfa::getOrCreateField(sim, dim, "ux0_N", nodes);
    tfa::getOrCreateField(sim, dim, "uy0_N", nodes);
    tfa::getOrCreateField(sim, dim, "uz0_N", nodes);
}

TF_Func bool RKiteration(tfa::Simulation &sim)
{
    CPUTimer_Beg("iteration");
    auto &M     = tfa::getMatrix(sim, "Id_NC");
    auto &GX    = tfa::getMatrix(sim, "GradX_CF");
    auto &GY    = tfa::getMatrix(sim, "GradY_CF");
    auto &GZ    = tfa::getMatrix(sim, "GradZ_CF");
    auto &ICF   = tfa::getMatrix(sim, "Interp_CF");
    auto &IFC   = tfa::getMatrix(sim, "Interp_FC");
    auto &ID_CN = tfa::getMatrix(sim, "Id_CN");

    auto &ux     = tfa::getField(sim, "ux_N");
    auto &uy     = tfa::getField(sim, "uy_N");
    auto &uz     = tfa::getField(sim, "uz_N");
    
    auto &ufx    = tfa::getField(sim, "ux_F");
    auto &ufy    = tfa::getField(sim, "uy_F");
    auto &ufz    = tfa::getField(sim, "uz_F");

    auto &upx    = tfa::getField(sim, "upredx_C");
    auto &upy    = tfa::getField(sim, "upredy_C");
    auto &upz    = tfa::getField(sim, "upredz_C");

    auto &msx    = tfa::getField(sim, "momSrcx_C");
    auto &msy    = tfa::getField(sim, "momSrcy_C");
    auto &msz    = tfa::getField(sim, "momSrcz_C");
    
    auto &p      = tfa::getField(sim, "P_C");
    static auto psolver = tfa::getSolver(sim, "Pressure_Solver");
    
    butcherTableau coefs = intScheme(sim.IOParamS["RKmethod"]);

    std::vector<double> b = coefs.b;
    std::vector<double> A = coefs.A;
    long unsigned int s = b.size();
 
    std::vector<tfa::Field*> convx(s);
    std::vector<tfa::Field*> diffx(s);
    std::vector<tfa::Field*> convy(s);
    std::vector<tfa::Field*> diffy(s);
    std::vector<tfa::Field*> convz(s);
    std::vector<tfa::Field*> diffz(s);  

    // Temporary variables.
    auto pSource = tfa::getTmpFieldS(sim, p);
    
    auto uxn = tfa::getTmpFieldS(sim,upx);
    auto uyn = tfa::getTmpFieldS(sim,upx);
    auto uzn = tfa::getTmpFieldS(sim,upx);

    tfa::oper_prod(M, ux, upx);
    tfa::oper_prod(M, uy, upy);
    tfa::oper_prod(M, uz, upz);

    tfa::oper_copy(upx,*uxn);
    tfa::oper_copy(upy,*uyn);
    tfa::oper_copy(upz,*uzn);

    double dt = sim.IOParamD["_TimeStep"];

    //SUBSTAGE 1
    
    // Predictor velocity should be computed here but for 1st stage it is not
    // required

    tfa::oper_prod(ICF, upx, ufx);
    tfa::oper_prod(ICF, upy, ufy);
    tfa::oper_prod(ICF, upz, ufz);
   
    auto &diffx0 = tfa::getOrCreateField(sim,"diffx_0",upx);
    auto &diffy0 = tfa::getOrCreateField(sim,"diffy_0",upx);
    auto &diffz0 = tfa::getOrCreateField(sim,"diffz_0",upx);

    auto &convx0 = tfa::getOrCreateField(sim,"convx_0",upx);
    auto &convy0 = tfa::getOrCreateField(sim,"convy_0",upx);
    auto &convz0 = tfa::getOrCreateField(sim,"convz_0",upx);

    tfa::oper_prod(diffx0,diffx0,0.0);
    tfa::oper_prod(diffy0,diffy0,0.0);
    tfa::oper_prod(diffz0,diffz0,0.0);
    tfa::oper_prod(convx0,convx0,0.0);
    tfa::oper_prod(convy0,convy0,0.0);
    tfa::oper_prod(convz0,convz0,0.0);

    diffx.at(0) = &diffx0;
    diffy.at(0) = &diffy0;
    diffz.at(0) = &diffz0;
    convx.at(0) = &convx0;
    convy.at(0) = &convy0;
    convz.at(0) = &convz0; 

    diffusive(ux,diffx0,sim);
    diffusive(uy,diffy0,sim);
    diffusive(uz,diffz0,sim);
    convective(ufx,convx0,sim);
    convective(ufy,convy0,sim);
    convective(ufz,convz0,sim);

    //predictor 2 stage

    for(long unsigned int i = 1; i<s; ++i){
      
      tfa::oper_copy(*uxn,upx);
      tfa::oper_copy(*uyn,upy);
      tfa::oper_copy(*uzn,upz);

      for(long unsigned int j=0; j<i; ++j){
        predictor(upx,upx,*(convx.at(j)),*(diffx.at(j)),msx,A.at(id(i+1,j+1)),dt);
        predictor(upy,upy,*(convy.at(j)),*(diffy.at(j)),msy,A.at(id(i+1,j+1)),dt);
        predictor(upz,upz,*(convz.at(j)),*(diffz.at(j)),msz,A.at(id(i+1,j+1)),dt);
      }
      
      tfa::oper_prod(ICF, upx, ufx);
      tfa::oper_prod(ICF, upy, ufy);
      tfa::oper_prod(ICF, upz, ufz);
    
      pRHS(*pSource,ufx,ufy,ufz,sim); 
      tfa::oper_solve(psolver, *pSource, p);
 
      projection(ufx,p,GX);
      projection(ufy,p,GY);
      projection(ufz,p,GZ);

      projection(upx,p,GX,IFC,sim);
      projection(upy,p,GY,IFC,sim);
      projection(upz,p,GZ,IFC,sim);
      
   // Map the velocity to the nodes.
      ID_CN_bocos(upx, ux, "ux_N", sim);
      ID_CN_bocos(upy, uy, "uy_N", sim);
      ID_CN_bocos(upz, uz, "uz_N", sim);

      //FINAL STAGE
 
      auto &convxi = tfa::getOrCreateField(sim,"convx_"+std::to_string(i),upx);
      auto &convyi = tfa::getOrCreateField(sim,"convy_"+std::to_string(i),upx);
      auto &convzi = tfa::getOrCreateField(sim,"convz_"+std::to_string(i),upx);
      auto &diffxi = tfa::getOrCreateField(sim,"diffx_"+std::to_string(i),upx);
      auto &diffyi = tfa::getOrCreateField(sim,"diffy_"+std::to_string(i),upx);
      auto &diffzi = tfa::getOrCreateField(sim,"diffz_"+std::to_string(i),upx);
 
      diffx.at(i) = &diffxi;
      diffy.at(i) = &diffyi;
      diffz.at(i) = &diffzi;
      convx.at(i) = &convxi;
      convy.at(i) = &convyi;
      convz.at(i) = &convzi;       

      tfa::oper_prod(diffxi,diffxi,0.0);
      tfa::oper_prod(diffyi,diffyi,0.0);
      tfa::oper_prod(diffzi,diffzi,0.0);
      tfa::oper_prod(convxi,convxi,0.0);
      tfa::oper_prod(convyi,convyi,0.0);
      tfa::oper_prod(convzi,convzi,0.0);

      diffusive(ux,diffxi,sim);
      diffusive(uy,diffyi,sim);
      diffusive(uz,diffzi,sim);
      convective(ufx,convxi,sim);
      convective(ufy,convyi,sim);
      convective(ufz,convzi,sim);
    }
    //printMaxVals("end of first stage. DIFFUSIVE FIELD",diffx,diffy,diffz);
    //printMaxVals("end of first stage. CONVECTIVE FIELD",convx,convy,convz);
    
    //Predictor velocity
    tfa::oper_copy(*uxn,upx);
    tfa::oper_copy(*uyn,upy);
    tfa::oper_copy(*uzn,upz);
    
    for(long unsigned int i = 0; i < s; ++i){ 
      predictor(upx,upx,*(convx.at(i)),*(diffx.at(i)),msx,b.at(i),dt);
      predictor(upy,upy,*(convy.at(i)),*(diffy.at(i)),msy,b.at(i),dt);
      predictor(upz,upz,*(convz.at(i)),*(diffz.at(i)),msz,b.at(i),dt);
    }

    tfa::oper_prod(ICF, upx, ufx);
    tfa::oper_prod(ICF, upy, ufy);
    tfa::oper_prod(ICF, upz, ufz);
    
    pRHS(*pSource,ufx,ufy,ufz,sim); 
    tfa::oper_solve(psolver, *pSource, p);
 
    projection(ufx,p,GX);
    projection(ufy,p,GY);
    projection(ufz,p,GZ);

    projection(upx,p,GX,IFC,sim);
    projection(upy,p,GY,IFC,sim);
    projection(upz,p,GZ,IFC,sim);
    
   // Map the velocity to the nodes.
    tfa::oper_prod(ID_CN, upx, ux);
    tfa::oper_prod(ID_CN, upy, uy);
    tfa::oper_prod(ID_CN, upz, uz);

    CPUTimer_End("iteration");
    return tfa::Iter_Continue;
}


TF_Func bool RKiteration_energy(tfa::Simulation &sim)
{
    CPUTimer_Beg("iteration");
    
    auto &M     = tfa::getMatrix(sim, "Id_NC");
    auto &GX    = tfa::getMatrix(sim, "GradX_CF");
    auto &GY    = tfa::getMatrix(sim, "GradY_CF");
    auto &GZ    = tfa::getMatrix(sim, "GradZ_CF");
    auto &ICF   = tfa::getMatrix(sim, "Interp_CF");
    auto &IFC   = tfa::getMatrix(sim, "Interp_FC");
    auto &ID_CN = tfa::getMatrix(sim, "Id_CN");
    auto &INF   = tfa::getMatrix(sim, "Interp_NF");
    
    auto &ux     = tfa::getField(sim, "ux_N");
    auto &uy     = tfa::getField(sim, "uy_N");
    auto &uz     = tfa::getField(sim, "uz_N");
    
    auto &T      = tfa::getField(sim, "T_N");
    
    auto &ufx    = tfa::getField(sim, "ux_F");
    auto &ufy    = tfa::getField(sim, "uy_F");
    auto &ufz    = tfa::getField(sim, "uz_F");

    auto &upx    = tfa::getField(sim, "upredx_C");
    auto &upy    = tfa::getField(sim, "upredy_C");
    auto &upz    = tfa::getField(sim, "upredz_C");

    auto &msx    = tfa::getField(sim, "momSrcx_C");
    auto &msy    = tfa::getField(sim, "momSrcy_C");
    auto &msz    = tfa::getField(sim, "momSrcz_C");
    auto msT    = tfa::getTmpFieldS(sim, msx);

    tfa::oper_setConst(*msT,0.0);
    
    auto &p      = tfa::getField(sim, "P_C");
    static auto psolver = tfa::getSolver(sim, "Pressure_Solver");
    butcherTableau coefs = intScheme(sim.IOParamS["RKmethod"]);

    double lambda = sim.IOParamD["lambda"]; //change when set properly

    std::vector<double> b = coefs.b;
    std::vector<double> A = coefs.A;
    long unsigned int s = b.size();
 
    std::vector<tfa::Field*> convx(s);
    std::vector<tfa::Field*> diffx(s);
    std::vector<tfa::Field*> convy(s);
    std::vector<tfa::Field*> diffy(s);
    std::vector<tfa::Field*> convz(s);
    std::vector<tfa::Field*> diffz(s);  

    std::vector<tfa::Field*> Trk(s); //stores the temperature at the cells for the Runge-Kutta substages
    
    // Temporary variables.
    auto pSource = tfa::getTmpFieldS(sim, p);
    
    auto uxn = tfa::getTmpFieldS(sim,upx);
    auto uyn = tfa::getTmpFieldS(sim,upx);
    auto uzn = tfa::getTmpFieldS(sim,upx);

    auto Tc = tfa::getTmpFieldS(sim,upx);
    auto Tf = tfa::getTmpFieldS(sim,ufx);
    
    tfa::oper_prod(M, ux, upx);
    tfa::oper_prod(M, uy, upy);
    tfa::oper_prod(M, uz, upz);

    tfa::oper_copy(upx,*uxn);
    tfa::oper_copy(upy,*uyn);
    tfa::oper_copy(upz,*uzn);

    double dt = sim.IOParamD["_TimeStep"];

    //SUBSTAGE 1
    
    // Predictor velocity should be computed here but for 1st stage it is not
    // required

    tfa::oper_prod(ICF, upx, ufx);
    tfa::oper_prod(ICF, upy, ufy);
    tfa::oper_prod(ICF, upz, ufz);
   
    auto &diffx0 = tfa::getOrCreateField(sim,"diffx_0",upx);
    auto &diffy0 = tfa::getOrCreateField(sim,"diffy_0",upx);
    auto &diffz0 = tfa::getOrCreateField(sim,"diffz_0",upx);

    auto &convx0 = tfa::getOrCreateField(sim,"convx_0",upx);
    auto &convy0 = tfa::getOrCreateField(sim,"convy_0",upx);
    auto &convz0 = tfa::getOrCreateField(sim,"convz_0",upx);

    auto &diffT = tfa::getOrCreateField(sim,"diffT",upx);
    auto &convT = tfa::getOrCreateField(sim,"convT",upx);
    auto &T0 = tfa::getOrCreateField(sim,"T_0",T);

    tfa::oper_prod(diffx0,diffx0,0.0);
    tfa::oper_prod(diffy0,diffy0,0.0);
    tfa::oper_prod(diffz0,diffz0,0.0);
    tfa::oper_prod(convx0,convx0,0.0);
    tfa::oper_prod(convy0,convy0,0.0);
    tfa::oper_prod(convz0,convz0,0.0);

    diffx.at(0) = &diffx0;
    diffy.at(0) = &diffy0;
    diffz.at(0) = &diffz0;
    convx.at(0) = &convx0;
    convy.at(0) = &convy0;
    convz.at(0) = &convz0; 

    diffusive(ux,diffx0,sim);
    diffusive(uy,diffy0,sim);
    diffusive(uz,diffz0,sim);
    convective(ufx,convx0,sim);
    convective(ufy,convy0,sim);
    convective(ufz,convz0,sim);

    tfa::oper_prod(T0,T0,0.0);
    Trk.at(0) = &T0;
    tfa::oper_copy(T,T0);

    tfa::oper_setConst(msx,0.0);
    tfa::oper_setConst(msz,0.0);
    
    if(sim.IOParamS["updateBoussinesq"] != "yes")
      generateSourceTerm(1.0,T,msy,sim);
    
    //predictor 2 stage

    for(long unsigned int i = 1; i<s; ++i){
      
      tfa::oper_copy(*uxn,upx);
      tfa::oper_copy(*uyn,upy);
      tfa::oper_copy(*uzn,upz);

      for(long unsigned int j=0; j<i; ++j){
        if(sim.IOParamS["updateBoussinesq"] == "yes")
          generateSourceTerm(1.0,*(Trk.at(j)),msy,sim); 
        
        predictor(upx,upx,*(convx.at(j)),*(diffx.at(j)),msx,A.at(id(i+1,j+1)),dt);
        predictor(upy,upy,*(convy.at(j)),*(diffy.at(j)),msy,A.at(id(i+1,j+1)),dt);
        predictor(upz,upz,*(convz.at(j)),*(diffz.at(j)),msz,A.at(id(i+1,j+1)),dt);
      }
      
      tfa::oper_prod(ICF, upx, ufx);
      tfa::oper_prod(ICF, upy, ufy);
      tfa::oper_prod(ICF, upz, ufz);
    
      pRHS(*pSource,ufx,ufy,ufz,sim); 
   
      CPUTimer_Beg("poisson");
      tfa::oper_solve(psolver, *pSource, p);
      CPUTimer_End("poisson");
      
      projection(ufx,p,GX);
      projection(ufy,p,GY);
      projection(ufz,p,GZ);

      projection(upx,p,GX,IFC,sim);
      projection(upy,p,GY,IFC,sim);
      projection(upz,p,GZ,IFC,sim);
      
   // Map the velocity to the nodes.
      ID_CN_bocos(upx, ux, "ux_N", sim);
      ID_CN_bocos(upy, uy, "uy_N", sim);
      ID_CN_bocos(upz, uz, "uz_N", sim);

      //ENERGY EQUATION
      
      tfa::oper_prod(M,T0,*Tc);
      
      for(long unsigned int j=0; j<i; ++j){
        tfa::oper_prod(INF,*(Trk.at(j)),*Tf);
        convective(*Tf,convT,sim);
        diffusive(*(Trk.at(j)),diffT,lambda,sim);
        predictor(*Tc,*Tc,convT,diffT,*msT,A.at(id(i+1,j+1)),dt);
      }
      
      auto &Ti = tfa::getOrCreateField(sim,"T_"+std::to_string(i),T0);

      Trk.at(i) = &Ti;
      tfa::oper_prod(Ti,Ti,0.0);
      ID_CN_bocos(*Tc,Ti,"T_N",sim);

      auto &convxi = tfa::getOrCreateField(sim,"convx_"+std::to_string(i),upx);
      auto &convyi = tfa::getOrCreateField(sim,"convy_"+std::to_string(i),upx);
      auto &convzi = tfa::getOrCreateField(sim,"convz_"+std::to_string(i),upx);
      auto &diffxi = tfa::getOrCreateField(sim,"diffx_"+std::to_string(i),upx);
      auto &diffyi = tfa::getOrCreateField(sim,"diffy_"+std::to_string(i),upx);
      auto &diffzi = tfa::getOrCreateField(sim,"diffz_"+std::to_string(i),upx);
 
      diffx.at(i) = &diffxi;
      diffy.at(i) = &diffyi;
      diffz.at(i) = &diffzi;
      convx.at(i) = &convxi;
      convy.at(i) = &convyi;
      convz.at(i) = &convzi;       

      tfa::oper_prod(diffxi,diffxi,0.0);
      tfa::oper_prod(diffyi,diffyi,0.0);
      tfa::oper_prod(diffzi,diffzi,0.0);
      tfa::oper_prod(convxi,convxi,0.0);
      tfa::oper_prod(convyi,convyi,0.0);
      tfa::oper_prod(convzi,convzi,0.0);

      diffusive(ux,diffxi,sim);
      diffusive(uy,diffyi,sim);
      diffusive(uz,diffzi,sim);
      convective(ufx,convxi,sim);
      convective(ufy,convyi,sim);
      convective(ufz,convzi,sim);
    }
    //printMaxVals("end of first stage. DIFFUSIVE FIELD",diffx,diffy,diffz);
    //printMaxVals("end of first stage. CONVECTIVE FIELD",convx,convy,convz);
    
    //Predictor velocity
    tfa::oper_copy(*uxn,upx);
    tfa::oper_copy(*uyn,upy);
    tfa::oper_copy(*uzn,upz);
    
    for(long unsigned int i = 0; i < s; ++i){ 
      if(sim.IOParamS["updateBoussinesq"] == "yes")
        generateSourceTerm(1.0,*(Trk.at(i)),msy,sim); 

      predictor(upx,upx,*(convx.at(i)),*(diffx.at(i)),msx,b.at(i),dt);
      predictor(upy,upy,*(convy.at(i)),*(diffy.at(i)),msy,b.at(i),dt);
      predictor(upz,upz,*(convz.at(i)),*(diffz.at(i)),msz,b.at(i),dt);
    }

    tfa::oper_prod(ICF, upx, ufx);
    tfa::oper_prod(ICF, upy, ufy);
    tfa::oper_prod(ICF, upz, ufz);
    
    pRHS(*pSource,ufx,ufy,ufz,sim);
    
    CPUTimer_Beg("poisson");
    tfa::oper_solve(psolver, *pSource, p);
    CPUTimer_End("poisson");

    projection(ufx,p,GX);
    projection(ufy,p,GY);
    projection(ufz,p,GZ);

    projection(upx,p,GX,IFC,sim);
    projection(upy,p,GY,IFC,sim);
    projection(upz,p,GZ,IFC,sim);
    
   // Map the velocity to the nodes.
    tfa::oper_prod(ID_CN, upx, ux);
    tfa::oper_prod(ID_CN, upy, uy);
    tfa::oper_prod(ID_CN, upz, uz);

    for(long unsigned int i = 0; i < s; ++i){ 
      tfa::oper_prod(INF,*(Trk.at(i)),*Tf);
      convective(*Tf,convT,sim);
      diffusive(*(Trk.at(i)),diffT,lambda,sim);
      predictor(*Tc,*Tc,convT,diffT,*msT,b.at(i),dt);
    }

    tfa::oper_prod(ID_CN, *Tc, T);

    CPUTimer_End("iteration");
    
    return tfa::Iter_Continue;
    //printMaxVals("end of time-step",ux,uy,uz);
}



