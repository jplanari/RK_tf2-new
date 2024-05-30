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
    TF_uAssert(tfa::hasField(sim, "ux"), "required field ux_N not defined");
    TF_uAssert(tfa::hasField(sim, "uy"), "required field uy_N not defined");
    TF_uAssert(tfa::hasField(sim, "uz"), "required field uz_N not defined");

    const auto &meshName = tfa::getMeshName(sim);
    const auto &msuf     = tfa::meshSuffix(meshName);
    const auto &faces    = "Faces"+msuf;
    const auto &nodes    = "Nodes"+msuf;
    
    const auto dim = tfa::getField(sim, "ux").dim;

    tfa::getOrCreateField(sim, dim, "momSrcx", nodes);
    tfa::getOrCreateField(sim, dim, "momSrcy", nodes);
    tfa::getOrCreateField(sim, dim, "momSrcz", nodes);

    tfa::getOrCreateField(sim, dim, "P",       nodes);

    tfa::getOrCreateField(sim, dim, "ux_F", faces);
    tfa::getOrCreateField(sim, dim, "uy_F", faces);
    tfa::getOrCreateField(sim, dim, "uz_F", faces);

    tfa::getOrCreateField(sim, dim, "ux0", nodes);
    tfa::getOrCreateField(sim, dim, "uy0", nodes);
    tfa::getOrCreateField(sim, dim, "uz0", nodes);
}

TF_Func bool RKiteration(tfa::Simulation &sim)
{
    CPUTimer_Beg("iteration");

    auto &GX    = tfa::getMatrix(sim, "GradX_NF");
    auto &GY    = tfa::getMatrix(sim, "GradY_NF");
    auto &GZ    = tfa::getMatrix(sim, "GradZ_NF");
    auto &IFN   = tfa::getMatrix(sim, "Interp_FN");
    auto &INF   = tfa::getMatrix(sim, "Interp_NF");
    
    auto &ux     = tfa::getField(sim, "ux");
    auto &uy     = tfa::getField(sim, "uy");
    auto &uz     = tfa::getField(sim, "uz");
    
    auto &ufx    = tfa::getField(sim, "ux_F");
    auto &ufy    = tfa::getField(sim, "uy_F");
    auto &ufz    = tfa::getField(sim, "uz_F");

    auto &msx    = tfa::getField(sim, "momSrcx");
    auto &msy    = tfa::getField(sim, "momSrcy");
    auto &msz    = tfa::getField(sim, "momSrcz");

    auto &p      = tfa::getField(sim, "P");
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
    
    auto uxn = tfa::getTmpFieldS(sim,ux);
    auto uyn = tfa::getTmpFieldS(sim,ux);
    auto uzn = tfa::getTmpFieldS(sim,ux);

    tfa::oper_copy(ux,*uxn);
    tfa::oper_copy(uy,*uyn);
    tfa::oper_copy(uz,*uzn);

    double dt = sim.IOParamD["_TimeStep"];

    //SUBSTAGE 1
    
    // Predictor velocity should be computed here but for 1st stage it is not
    // required

    tfa::oper_prod(INF, ux, ufx);
    tfa::oper_prod(INF, uy, ufy);
    tfa::oper_prod(INF, uz, ufz);
   
    auto &diffx0 = tfa::getOrCreateField(sim,"diffx_0",ux);
    auto &diffy0 = tfa::getOrCreateField(sim,"diffy_0",ux);
    auto &diffz0 = tfa::getOrCreateField(sim,"diffz_0",ux);

    auto &convx0 = tfa::getOrCreateField(sim,"convx_0",ux);
    auto &convy0 = tfa::getOrCreateField(sim,"convy_0",ux);
    auto &convz0 = tfa::getOrCreateField(sim,"convz_0",ux);

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
     
      tfa::oper_copy(*uxn,ux);
      tfa::oper_copy(*uyn,uy);
      tfa::oper_copy(*uzn,uz);

      for(long unsigned int j=0; j<i; ++j){
        predictor(ux,ux,*(convx.at(j)),*(diffx.at(j)),msx,A.at(id(i+1,j+1)),dt);
        predictor(uy,uy,*(convy.at(j)),*(diffy.at(j)),msy,A.at(id(i+1,j+1)),dt);
        predictor(uz,uz,*(convz.at(j)),*(diffz.at(j)),msz,A.at(id(i+1,j+1)),dt);
      }
      
      tfa::oper_prod(INF, ux, ufx);
      tfa::oper_prod(INF, uy, ufy);
      tfa::oper_prod(INF, uz, ufz);
    
      pRHS(*pSource,ufx,ufy,ufz,sim); 
   
      CPUTimer_Beg("poisson");
      tfa::oper_solve(psolver, *pSource, p);
      CPUTimer_End("poisson");
      
      projection(ufx,p,GX);
      projection(ufy,p,GY);
      projection(ufz,p,GZ);

      projection(ux,p,GX,IFN,sim);
      projection(uy,p,GY,IFN,sim);
      projection(uz,p,GZ,IFN,sim);
      
   // Map the velocity to the nodes.
      innerRK_bocos(ux,"ux", sim);
      innerRK_bocos(uy,"uy", sim);
      innerRK_bocos(uz,"uz", sim);

      //FINAL STAGE
 
      auto &convxi = tfa::getOrCreateField(sim,"convx_"+std::to_string(i),ux);
      auto &convyi = tfa::getOrCreateField(sim,"convy_"+std::to_string(i),ux);
      auto &convzi = tfa::getOrCreateField(sim,"convz_"+std::to_string(i),ux);
      auto &diffxi = tfa::getOrCreateField(sim,"diffx_"+std::to_string(i),ux);
      auto &diffyi = tfa::getOrCreateField(sim,"diffy_"+std::to_string(i),ux);
      auto &diffzi = tfa::getOrCreateField(sim,"diffz_"+std::to_string(i),ux);
 
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
    
    tfa::oper_copy(*uxn,ux);
    tfa::oper_copy(*uyn,uy);
    tfa::oper_copy(*uzn,uz);
    
    for(long unsigned int i = 0; i < s; ++i){ 
      predictor(ux,ux,*(convx.at(i)),*(diffx.at(i)),msx,b.at(i),dt);
      predictor(uy,uy,*(convy.at(i)),*(diffy.at(i)),msy,b.at(i),dt);
      predictor(uz,uz,*(convz.at(i)),*(diffz.at(i)),msz,b.at(i),dt);
    }

    tfa::oper_prod(INF, ux, ufx);
    tfa::oper_prod(INF, uy, ufy);
    tfa::oper_prod(INF, uz, ufz);
    
    pRHS(*pSource,ufx,ufy,ufz,sim);
    
    CPUTimer_Beg("poisson");
    tfa::oper_solve(psolver, *pSource, p);
    CPUTimer_End("poisson");

    projection(ufx,p,GX);
    projection(ufy,p,GY);
    projection(ufz,p,GZ);

    projection(ux,p,GX,IFN,sim);
    projection(uy,p,GY,IFN,sim);
    projection(uz,p,GZ,IFN,sim);
   
    CPUTimer_End("iteration");
    return tfa::Iter_Continue;
}


TF_Func bool RKiteration_energy(tfa::Simulation &sim)
{
    CPUTimer_Beg("iteration");
    
    tfa::info("I am in");
    auto &GX    = tfa::getMatrix(sim, "GradX_NF");
    auto &GY    = tfa::getMatrix(sim, "GradY_NF");
    auto &GZ    = tfa::getMatrix(sim, "GradZ_NF");
    auto &IFN   = tfa::getMatrix(sim, "Interp_FN");
    auto &INF   = tfa::getMatrix(sim, "Interp_NF");
    
    auto &ux     = tfa::getField(sim, "ux");
    auto &uy     = tfa::getField(sim, "uy");
    auto &uz     = tfa::getField(sim, "uz");
    
    auto &T      = tfa::getField(sim, "T");
    
    auto &ufx    = tfa::getField(sim, "ux_F");
    auto &ufy    = tfa::getField(sim, "uy_F");
    auto &ufz    = tfa::getField(sim, "uz_F");

    auto &msx    = tfa::getField(sim, "momSrcx");
    auto &msy    = tfa::getField(sim, "momSrcy");
    auto &msz    = tfa::getField(sim, "momSrcz");
    auto msT    = tfa::getTmpFieldS(sim, msx);

    tfa::oper_setConst(*msT,0.0);
    
    
    auto &p      = tfa::getField(sim, "P");
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
    
    auto uxn = tfa::getTmpFieldS(sim,ux);
    auto uyn = tfa::getTmpFieldS(sim,ux);
    auto uzn = tfa::getTmpFieldS(sim,ux);
    
    auto Tf = tfa::getTmpFieldS(sim,ufx);

    tfa::oper_copy(ux,*uxn);
    tfa::oper_copy(uy,*uyn);
    tfa::oper_copy(uz,*uzn);

    double dt = sim.IOParamD["_TimeStep"];
    

    //SUBSTAGE 1
    
    // Predictor velocity should be computed here but for 1st stage it is not
    // required

    tfa::oper_prod(INF, ux, ufx);
    tfa::oper_prod(INF, uy, ufy);
    tfa::oper_prod(INF, uz, ufz);
   
    auto &diffx0 = tfa::getOrCreateField(sim,"diffx_0",ux);
    auto &diffy0 = tfa::getOrCreateField(sim,"diffy_0",ux);
    auto &diffz0 = tfa::getOrCreateField(sim,"diffz_0",ux);

    auto &convx0 = tfa::getOrCreateField(sim,"convx_0",ux);
    auto &convy0 = tfa::getOrCreateField(sim,"convy_0",ux);
    auto &convz0 = tfa::getOrCreateField(sim,"convz_0",ux);

    auto &diffT = tfa::getOrCreateField(sim,"diffT",T);
    auto &convT = tfa::getOrCreateField(sim,"convT",T);
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
      
      tfa::oper_copy(*uxn,ux);
      tfa::oper_copy(*uyn,uy);
      tfa::oper_copy(*uzn,uz);

      for(long unsigned int j=0; j<i; ++j){
        
        if(sim.IOParamS["updateBoussinesq"] == "yes")
          generateSourceTerm(1.0,*(Trk.at(j)),msy,sim); 
        
        predictor(ux,ux,*(convx.at(j)),*(diffx.at(j)),msx,A.at(id(i+1,j+1)),dt);
        predictor(uy,uy,*(convy.at(j)),*(diffy.at(j)),msy,A.at(id(i+1,j+1)),dt);
        predictor(uz,uz,*(convz.at(j)),*(diffz.at(j)),msz,A.at(id(i+1,j+1)),dt);
      
      }
      
      tfa::oper_prod(INF, ux, ufx);
      tfa::oper_prod(INF, uy, ufy);
      tfa::oper_prod(INF, uz, ufz);
    
      pRHS(*pSource,ufx,ufy,ufz,sim); 
   
      CPUTimer_Beg("poisson");
      tfa::oper_solve(psolver, *pSource, p);
      CPUTimer_End("poisson");
      
      projection(ufx,p,GX);
      projection(ufy,p,GY);
      projection(ufz,p,GZ);

      projection(ux,p,GX,IFN,sim);
      projection(uy,p,GY,IFN,sim);
      projection(uz,p,GZ,IFN,sim);
      
   // Map the velocity to the nodes.
      innerRK_bocos(ux,"ux", sim);
      innerRK_bocos(uy,"uy", sim);
      innerRK_bocos(uz,"uz", sim);

      //ENERGY EQUATION
      
      auto &Ti = tfa::getOrCreateField(sim,"T_"+std::to_string(i),T0);
      Trk.at(i) = &Ti;
      
      tfa::oper_copy(T0,Ti);

      for(long unsigned int j=0; j<i; ++j){
        tfa::oper_prod(INF,*(Trk.at(j)),*Tf);
        convective(*Tf,convT,sim);
        diffusive(*(Trk.at(j)),diffT,lambda,sim);
        predictor(Ti,Ti,convT,diffT,*msT,A.at(id(i+1,j+1)),dt);
      }
      innerRK_bocos(Ti,"T",sim);

      auto &convxi = tfa::getOrCreateField(sim,"convx_"+std::to_string(i),ux);
      auto &convyi = tfa::getOrCreateField(sim,"convy_"+std::to_string(i),ux);
      auto &convzi = tfa::getOrCreateField(sim,"convz_"+std::to_string(i),ux);
      auto &diffxi = tfa::getOrCreateField(sim,"diffx_"+std::to_string(i),ux);
      auto &diffyi = tfa::getOrCreateField(sim,"diffy_"+std::to_string(i),ux);
      auto &diffzi = tfa::getOrCreateField(sim,"diffz_"+std::to_string(i),ux);
 
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
    tfa::oper_copy(*uxn,ux);
    tfa::oper_copy(*uyn,uy);
    tfa::oper_copy(*uzn,uz);
    
    for(long unsigned int i = 0; i < s; ++i){ 
      if(sim.IOParamS["updateBoussinesq"] == "yes")
        generateSourceTerm(1.0,*(Trk.at(i)),msy,sim); 

      predictor(ux,ux,*(convx.at(i)),*(diffx.at(i)),msx,b.at(i),dt);
      predictor(uy,uy,*(convy.at(i)),*(diffy.at(i)),msy,b.at(i),dt);
      predictor(uz,uz,*(convz.at(i)),*(diffz.at(i)),msz,b.at(i),dt);
    }

    tfa::oper_prod(INF, ux, ufx);
    tfa::oper_prod(INF, uy, ufy);
    tfa::oper_prod(INF, uz, ufz);
    
    pRHS(*pSource,ufx,ufy,ufz,sim);
    
    CPUTimer_Beg("poisson");
    tfa::oper_solve(psolver, *pSource, p);
    CPUTimer_End("poisson");

    projection(ufx,p,GX);
    projection(ufy,p,GY);
    projection(ufz,p,GZ);

    projection(ux,p,GX,IFN,sim);
    projection(uy,p,GY,IFN,sim);
    projection(uz,p,GZ,IFN,sim);
    
    tfa::oper_copy(T0,T);  
  
    for(long unsigned int i = 0; i < s; ++i){ 
      tfa::oper_prod(INF,*(Trk.at(i)),*Tf);
      convective(*Tf,convT,sim);
      diffusive(*(Trk.at(i)),diffT,lambda,sim);
      predictor(T,T,convT,diffT,*msT,b.at(i),dt);
    }

    CPUTimer_End("iteration");
    
    return tfa::Iter_Continue;
    //printMaxVals("end of time-step",ux,uy,uz);
}



