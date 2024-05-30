tfa::Field &calcFaceStaggeredVolumes(tfa::Simulation &sim)
{
    if(tfa::hasField(sim, "_FaceVolume")) return tfa::getField(sim, "_FaceVolume");
    
    auto dim = tfa::getField(sim,"ux").dim;
    auto &result = tfa::newField(sim, dim, "_FaceVolume", "Faces");
    
    uint32_t resultSize = tfa::getNumEntries(result);
    std::vector<double> buffer(resultSize);

    if (tfa::hasSMesh(sim))
    {
      const auto &m = tfa::getSMesh(sim);
      for (uint32_t it = 0; it < resultSize/dim; it++)
      {
        const auto &F = tfa::getFaceIJKFromId(m,it);
        for(uint32_t d=0; d<dim; d++)
        buffer[it*dim+d] = tfa::calcFaceStaggeredVol(m,F);
      }
    }

    tfa::oper_setData(result, buffer.data());
    return result;
}

tfa::Field &faceCellDist(tfa::Simulation &sim)
{
  if(tfa::hasField(sim,"_FaceCellDist")) return tfa::getField(sim,"_FaceCellDist");

  auto dim = tfa::getField(sim,"ux").dim;
  auto &result = tfa::newField(sim, dim, "_FaceCellDist", "Faces");
  uint32_t resultSize = tfa::getNumEntries(result);
  std::vector<double> buffer(resultSize);

  if(tfa::hasSMesh(sim))
  {
    const auto &m = tfa::getSMesh(sim);
    for (uint32_t it = 0; it < resultSize/dim; it++)
    {
      const auto F = tfa::getFaceIJKFromId(m, it);
      const auto n = tfa::calcFaceNormal(m, F);
      if (tfa::isBoundaryFace(m, F))
      {
        const auto A = tfa::calcCellDim(m, F);
        for(uint32_t d=0; d<dim; d++)
            buffer[it*dim+d] = 0.5*fabs(A*n);
      }
      else
      {
        const auto O = tfa::getOtherCellIJKFromFace(m, F);
        for(uint32_t d=0; d<dim; d++)
            buffer[it*dim+d] = tfa::calcCellDist(m, F.i, F.j, F.k, O.i, O.j, O.k);
      }
    } 
  }

  tfa::oper_setData(result, buffer.data());
  return result;
}



tfa::IJK returnNbIJK(tfa::IJK CNB, uint32_t f)
{
  tfa::IJK C = CNB;
  switch(f){
    case 0:
      C.i--;
      break;
    case 1:
      C.i++;
      break;
    case 2:
      C.j--;
      break;
    case 3:
      C.j++;
      break;
    case 4:
      C.k--;
      break;
    case 5:
      C.k++;
      break;
  }
  return C;
}

std::vector<tfa::IJK> getFaceNbCells(tfa::IJK F, const tfa::SMesh &m)
{
  std::vector<tfa::IJK> result;
  result.push_back(F);
  auto f = F.f;
  //If face F is not a boundary face:
  if(!tfa::isBoundaryFace(m,F)){
     result.push_back(returnNbIJK(F,(uint32_t)f));}
  else
  {
    if(f%2!=0){ f--; result.at(0) = returnNbIJK(F,(uint32_t)f);}
    else result.at(0) = F;
  }

  return result; 
}
