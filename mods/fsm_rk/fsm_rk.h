void predictor(tfa::Field &up, tfa::Field &un, tfa::Field &conv, tfa::Field &diff, tfa::Field &ms, double coef, double dt)
{
    tfa::oper_copy(un, up);
    tfa::oper_axpy(diff, up, dt*coef, 1.0);
    tfa::oper_axpy(conv, up, -dt*coef, 1.0);
    tfa::oper_axpy(ms, up, dt*coef, 1.0);
}

void projection(tfa::Field &uf, tfa::Field &p, tfa::Matrix &G)
{
    tfa::oper_prod(G,p,uf,-1.0,1.0);
}

void projection(tfa::Field &up, tfa::Field &p, tfa::Matrix &G, tfa::Matrix &IFC, tfa::Simulation &sim)
{
    auto Gp = tfa::getTmpFieldS(sim,up.dim,"Faces");
    tfa::oper_prod(G,p,*Gp);
    tfa::oper_prod(IFC,*Gp,up,-1.0,1.0);
}

void pRHS(tfa::Field &pS, tfa::Field &ux, tfa::Field &uy, tfa::Field &uz, tfa::Simulation &sim)
{
    auto &Omega = tfa::getField(sim, "Omega_C");
    auto &DX    = tfa::getMatrix(sim, "DivX_FC");
    auto &DY    = tfa::getMatrix(sim, "DivY_FC");
    auto &DZ    = tfa::getMatrix(sim, "DivZ_FC");
    
    tfa::oper_prod(DX, ux, pS);
    tfa::oper_prod(DY, uy, pS, 1.0, 1.0);
    tfa::oper_prod(DZ, uz, pS, 1.0, 1.0);
    tfa::oper_prod(Omega, pS, -1.0);
}

void ID_CN_bocos(tfa::Field &phiC, tfa::Field &phiN, std::string fieldName, tfa::Simulation &sim)
{
  std::string bndsrcName = fieldName + "_BndSrc";
  std::string bndapplyName = fieldName + "_BndApply_NN";

  auto aux = tfa::getTmpFieldS(sim, phiN);
  auto &phi_bnd = tfa::getField(sim, bndsrcName);
  auto &MB = tfa::getMatrix(sim, bndapplyName);
  auto &cm = tfa::meshNodeCellMask(sim, phiN.dim);
  auto &ID_CN = tfa::getMatrix(sim, "Id_CN");
  
  tfa::oper_prod(ID_CN, phiC, *aux);
  tfa::oper_copy(phi_bnd, phiN);
  tfa::oper_fmadd(cm, *aux, phiN, 1.0, 1.0);
  tfa::oper_prod(MB, *aux, phiN, 1.0, 1.0);
}
