void generateSourceTerm(double alpha, tfa::Field &T, tfa::Field &source, tfa::Simulation &sim)
{
  tfa::oper_axpy(T,source,alpha,0.0);
}


