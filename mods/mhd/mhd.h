void generateSourceTerm(tfa::Field &Jx, tfa::Field &Jy, tfa::Field &Jz:w, tfa::Field &source, tfa::Simulation &sim)
{
  auto &Tc = tfa::getTmpFieldS(sim,source);
  auto &M = tfa::getMatrix(sim, "Id_NC");
  tfa::oper_prod(M,T,Tc);
  tfa::oper_axpy(Tc,source,alpha,0.0);
}


