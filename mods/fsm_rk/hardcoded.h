#define id(eps,sig) (eps-1)*(eps-2)/2+(sig-1)

TF_Func bool EulerIteration(tfa::Simulation &sim)
{
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
   
    // Temporary variables.
    auto &pSource = tfa::getTmpFieldS(sim, p);
    
    auto &uxn = tfa::getTmpFieldS(sim,upx);
    auto &uyn = tfa::getTmpFieldS(sim,upx);
    auto &uzn = tfa::getTmpFieldS(sim,upx);

    tfa::oper_prod(M, ux, upx);
    tfa::oper_prod(M, uy, upy);
    tfa::oper_prod(M, uz, upz);

    tfa::oper_copy(upx,uxn);
    tfa::oper_copy(upy,uyn);
    tfa::oper_copy(upz,uzn);

    double dt = sim.IOParamD["_TimeStep"];

    //SUBSTAGE 1
    
    // Predictor velocity should be computed here but for 1st stage it is not
    // required

    //printMaxVals("end of first stage",ux,uy,uz);

    //FINAL STAGE
    
    auto &diffx = diffusive(ux,0,sim);
    auto &diffy = diffusive(uy,1,sim);
    auto &diffz = diffusive(uz,2,sim);
    auto &convx = convective(ufx,0,sim);
    auto &convy = convective(ufy,1,sim);
    auto &convz = convective(ufz,2,sim);



    //printMaxVals("end of first stage. DIFFUSIVE FIELD",diffx,diffy,diffz);
    //printMaxVals("end of first stage. CONVECTIVE FIELD",convx,convy,convz);
    
    //Predictor velocity
    
    upx = predictor(uxn,convx,diffx,msx,dt,1.0,sim);
    upy = predictor(uyn,convy,diffy,msy,dt,1.0,sim);
    upz = predictor(uzn,convz,diffz,msz,dt,1.0,sim);
   
    tfa::oper_prod(ICF, upx, ufx);
    tfa::oper_prod(ICF, upy, ufy);
    tfa::oper_prod(ICF, upz, ufz);
    
    pSource = pRHS(ufx,ufy,ufz,sim); 
    tfa::oper_solve(psolver, pSource, p);
 
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

    //printMaxVals("end of time-step",ux,uy,uz);

    return tfa::Iter_Continue;
}

TF_Func bool RK2Iteration(tfa::Simulation &sim)
{
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

   
    // Temporary variables.
    auto &pSource = tfa::getTmpFieldS(sim, p);
    
    auto &uxn = tfa::getTmpFieldS(sim,upx);
    auto &uyn = tfa::getTmpFieldS(sim,upx);
    auto &uzn = tfa::getTmpFieldS(sim,upx);

    tfa::oper_prod(M, ux, upx);
    tfa::oper_prod(M, uy, upy);
    tfa::oper_prod(M, uz, upz);

    tfa::oper_copy(upx,uxn);
    tfa::oper_copy(upy,uyn);
    tfa::oper_copy(upz,uzn);

    double dt = sim.IOParamD["_TimeStep"];

    //SUBSTAGE 1
    
    // Predictor velocity should be computed here but for 1st stage it is not
    // required

    tfa::oper_prod(ICF, upx, ufx);
    tfa::oper_prod(ICF, upy, ufy);
    tfa::oper_prod(ICF, upz, ufz);
    
    //SUBSTAGE 2
    auto &diffx  = diffusive(ux,0,sim);
    auto &diffy  = diffusive(uy,1,sim);
    auto &diffz  = diffusive(uz,2,sim);
    auto &convx  = convective(ufx,0,sim);
    auto &convy  = convective(ufy,1,sim);
    auto &convz  = convective(ufz,2,sim);

    //predictor 2 stage

    upx = predictor(uxn,convx,diffx,msx,1.0,dt,sim);
    upy = predictor(uyn,convy,diffy,msy,1.0,dt,sim);
    upz = predictor(uzn,convz,diffz,msz,1.0,dt,sim);
   
    tfa::oper_prod(ICF, upx, ufx);
    tfa::oper_prod(ICF, upy, ufy);
    tfa::oper_prod(ICF, upz, ufz);
    
    pSource = pRHS(ufx,ufy,ufz,sim); 
    tfa::oper_solve(psolver, pSource, p);
 
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

    //FINAL STAGE
 
    auto &diffx2  = diffusive(ux,0,sim);
    auto &diffy2  = diffusive(uy,1,sim);
    auto &diffz2  = diffusive(uz,2,sim);
    auto &convx2  = convective(ufx,0,sim);
    auto &convy2  = convective(ufy,1,sim);
    auto &convz2  = convective(ufz,2,sim);
    
    //printMaxVals("end of first stage. DIFFUSIVE FIELD",diffx,diffy,diffz);
    //printMaxVals("end of first stage. CONVECTIVE FIELD",convx,convy,convz);
    
    //Predictor velocity
    
    upx = predictor(uxn,convx,diffx,msx,0.5,dt,sim);
    upy = predictor(uyn,convy,diffy,msy,0.5,dt,sim);
    upz = predictor(uzn,convz,diffz,msz,0.5,dt,sim);

    upx = predictor(upx,convx2,diffx2,msx,0.5,dt,sim);
    upy = predictor(upy,convy2,diffy2,msy,0.5,dt,sim);
    upz = predictor(upz,convz2,diffz2,msz,0.5,dt,sim);

    tfa::oper_prod(ICF, upx, ufx);
    tfa::oper_prod(ICF, upy, ufy);
    tfa::oper_prod(ICF, upz, ufz);
    
    pSource = pRHS(ufx,ufy,ufz,sim); 
    tfa::oper_solve(psolver, pSource, p);
 
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

    //printMaxVals("end of time-step",ux,uy,uz);

    return tfa::Iter_Continue;
}

TF_Func bool RK2Iteration_vector(tfa::Simulation &sim)
{
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
    int s = b.size();
 
    std::vector<tfa::Field*> convx(s);
    std::vector<tfa::Field*> diffx(s);
    std::vector<tfa::Field*> convy(s);
    std::vector<tfa::Field*> diffy(s);
    std::vector<tfa::Field*> convz(s);
    std::vector<tfa::Field*> diffz(s);  

    // Temporary variables.
    auto &pSource = tfa::getTmpFieldS(sim, p);
    
    auto &uxn = tfa::getTmpFieldS(sim,upx);
    auto &uyn = tfa::getTmpFieldS(sim,upx);
    auto &uzn = tfa::getTmpFieldS(sim,upx);

    tfa::oper_prod(M, ux, upx);
    tfa::oper_prod(M, uy, upy);
    tfa::oper_prod(M, uz, upz);

    tfa::oper_copy(upx,uxn);
    tfa::oper_copy(upy,uyn);
    tfa::oper_copy(upz,uzn);

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

    diffx0  = diffusive(ux,0,sim);
    diffy0  = diffusive(uy,1,sim);
    diffz0  = diffusive(uz,2,sim);
    convx0  = convective(ufx,0,sim);
    convy0  = convective(ufy,1,sim);
    convz0  = convective(ufz,2,sim);

    //predictor 2 stage

    for(int i = 1; i<s; ++i){
      
      tfa::oper_copy(uxn,upx);
      tfa::oper_copy(uyn,upy);
      tfa::oper_copy(uzn,upz);

      for(int j=0; j<i; ++j){
        upx = predictor(upx,*(convx.at(j)),*(diffx.at(j)),msx,A.at(id(i+1,j+1)),dt,sim);
        upy = predictor(upy,*(convy.at(j)),*(diffy.at(j)),msy,A.at(id(i+1,j+1)),dt,sim);
        upz = predictor(upz,*(convz.at(j)),*(diffz.at(j)),msz,A.at(id(i+1,j+1)),dt,sim);
      }
      
      tfa::oper_prod(ICF, upx, ufx);
      tfa::oper_prod(ICF, upy, ufy);
      tfa::oper_prod(ICF, upz, ufz);
    
      pSource = pRHS(ufx,ufy,ufz,sim); 
      tfa::oper_solve(psolver, pSource, p);
 
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

      diffxi  = diffusive(ux,0,sim);
      diffyi  = diffusive(uy,1,sim);
      diffzi  = diffusive(uz,2,sim);
      convxi  = convective(ufx,0,sim);
      convyi  = convective(ufy,1,sim);
      convzi  = convective(ufz,2,sim);
    }
    //printMaxVals("end of first stage. DIFFUSIVE FIELD",diffx,diffy,diffz);
    //printMaxVals("end of first stage. CONVECTIVE FIELD",convx,convy,convz);
    
    //Predictor velocity
    tfa::oper_copy(uxn,upx);
    tfa::oper_copy(uyn,upy);
    tfa::oper_copy(uzn,upz);
    
    for(int i = 0; i < s; ++i){ 
      upx = predictor(upx,*(convx.at(i)),*(diffx.at(i)),msx,b.at(i),dt,sim);
      upy = predictor(upy,*(convy.at(i)),*(diffy.at(i)),msy,b.at(i),dt,sim);
      upz = predictor(upz,*(convz.at(i)),*(diffz.at(i)),msz,b.at(i),dt,sim);
    }

    tfa::oper_prod(ICF, upx, ufx);
    tfa::oper_prod(ICF, upy, ufy);
    tfa::oper_prod(ICF, upz, ufz);
    
    pSource = pRHS(ufx,ufy,ufz,sim); 
    tfa::oper_solve(psolver, pSource, p);
 
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

    //printMaxVals("end of time-step",ux,uy,uz);

    return tfa::Iter_Continue;
}
