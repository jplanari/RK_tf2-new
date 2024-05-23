std::vector<std::string> multiRK_method(tfa::Simulation &sim)
{
  uint64_t Nmax = 9;
  std::vector<std::string> pool;
  std::string fcn_name, base_name = "RKmethod_";
  for(uint64_t i = 1; i <= Nmax; ++i)
  {
    fcn_name = sim.IOParamS[base_name + std::to_string(i)];
    if(!fcn_name.empty())
      pool.push_back(fcn_name);
    else
      break;
  }
  return pool;
}
