#include "meshFcs.h"

double computeSFaceContributionGershgorin_real(tfa::SMesh &m, tfa::IJK F, tfa::IJK C, tfa::IJK CNB)
{
    auto cntr = tfa::calcCellCentroid(m,C);
    auto cntr_nb = tfa::calcCellCentroid(m,CNB);
    tfa::V3 dvec = cntr-cntr_nb;
    auto fnv = tfa::calcFaceNormal(m,F);
    double area = fnv*tfa::calcCellFacesArea(m,C);

    return fabs(area/(dvec*fnv));
}

double computeFaceContributionGershgorin_imag(tfa::V3 fnv, tfa::V3 uf, double Af)
{
    return 0.5*fabs(fnv*uf)*Af;
}

void computeImagEV_Gershgorin(tfa::Simulation &sim)
{
    auto &ux_N = tfa::getField(sim,"ux_N");
    auto &uy_N = tfa::getField(sim,"uy_N");
    auto &uz_N = tfa::getField(sim,"uz_N");

    auto &INF = tfa::getMatrix(sim,"Interp_NF");
    
    auto uxf = tfa::getTmpFieldS(sim,ux_N.dim,"Faces");
    auto uyf = tfa::getTmpFieldS(sim,ux_N.dim,"Faces");
    auto uzf = tfa::getTmpFieldS(sim,ux_N.dim,"Faces");

    tfa::oper_prod(INF,ux_N,*uxf);
    tfa::oper_prod(INF,uy_N,*uyf);
    tfa::oper_prod(INF,uz_N,*uzf);

    double u,v,w;
    tfa::V3 uf;
    double gershContImag;
    double evi=0.0;
	    tfa::SMesh &m = tfa::getSMesh(sim);
      tfa::IJK F,C;
	    uint32_t offset = (uint32_t)m.topo.firstOwnedFaces[tfa::mpiRank()];
      uint32_t f;
	    tfa::V3 cfa,fnv;
	    for(uint32_t c = 0; c < tfa::getNumCells(m); ++c){
		
		    gershContImag = 0.0;
		    C = tfa::getCellIJKFromId(m,c);
		    F = C;
		
	      cfa = tfa::calcCellFacesArea(m,C);
		
		    for(int fc = 0; fc < 6; fc++){
			  
          if(hasNBCell(m,C,(uint32_t)fc)){
            F.f = fc;
				    f = (uint32_t)tfa::getFaceIndex(m,C.i,C.j,C.k,fc)-offset;
				    fnv = tfa::calcFaceNormal(m,F);
			
				    u = (*uxf).array[TO_LOC(f)];
				    v = (*uyf).array[TO_LOC(f)];
				    w = (*uzf).array[TO_LOC(f)];

				    uf = {u,v,w};
				    gershContImag += computeFaceContributionGershgorin_imag(fnv,uf,fnv*cfa)/tfa::calcCellVolume(m,C);			
			    }
		    }
		    evi = std::max(evi,gershContImag);
	    }
    sim.IOParamD["_EVimag"] = tfa::allReduce(evi, MPI_MAX);
}

TF_Func void computeEV_Gershgorin(tfa::Simulation &sim)
{
    auto &ux_N = tfa::getField(sim,"ux_N");
    auto &uy_N = tfa::getField(sim,"uy_N");
    auto &uz_N = tfa::getField(sim,"uz_N");
    auto &INF = tfa::getMatrix(sim,"Interp_NF");
    
    auto ux_f = tfa::getTmpFieldS(sim,ux_N.dim,"Faces");
    auto uy_f = tfa::getTmpFieldS(sim,ux_N.dim,"Faces");
    auto uz_f = tfa::getTmpFieldS(sim,ux_N.dim,"Faces");

    tfa::oper_prod(INF,ux_N,*ux_f);
    tfa::oper_prod(INF,uy_N,*uy_f);
    tfa::oper_prod(INF,uz_N,*uz_f);

    double kinVisc = sim.IOParamD["kinVisc"];

    double gershContImag;
    double gershContReal;

    double evr=0.0,evi=0.0;

    double u,v,w;
    tfa::V3 uf;

	    tfa::SMesh &m = tfa::getSMesh(sim);
      tfa::IJK F,C,CNB;
	    uint32_t offset = (uint32_t)m.topo.firstOwnedFaces[tfa::mpiRank()];
      uint32_t f;
	    tfa::V3 fnv, cfa;
	    for(uint32_t c = 0; c < tfa::getNumCells(m); ++c) {
		    gershContImag = 0.0;
		    gershContReal = 0.0;
		
		    C = tfa::getCellIJKFromId(m,c);
		    F = C;
		    CNB = tfa::getCellIJKFromId(m,c);
		    cfa = tfa::calcCellFacesArea(m,C);
		    for(int fc = 0; fc < 6; fc++){
			    if(hasNBCell(m,C,(uint32_t)fc)){
				    f = (uint32_t)tfa::getFaceIndex(m,C.i,C.j,C.k,fc)-offset;
				    F.f = fc;
				    CNB = C;
				    fnv = tfa::calcFaceNormal(m,F);
				
				    switch(fc){
				      case 0:
				    	  CNB.i = C.i-1;
					      break;
				      case 1:
				    	  CNB.i = C.i+1;
					      break;
				      case 2:
				        CNB.j = C.j-1;
					      break;
				      case 3:
				        CNB.j = C.j+1;
					      break;
				      case 4:
				        CNB.k = C.k-1;
					      break;
				      case 5:
					      CNB.k = C.k+1;
					      break;
				    }
			
				    u = (*ux_f).array[TO_LOC(f)];
				    v = (*uy_f).array[TO_LOC(f)];
				    w = (*uz_f).array[TO_LOC(f)];

				    uf = {u,v,w};
				    gershContReal += 2.0*kinVisc*computeSFaceContributionGershgorin_real(m,F,C,CNB)/tfa::calcCellVolume(m,C);
				    gershContImag += computeFaceContributionGershgorin_imag(fnv,uf,fnv*cfa)/tfa::calcCellVolume(m,C);			
			    }
		    }
        evr = std::max(evr,gershContReal);
		    evi = std::max(evi,gershContImag);
	    }
    sim.IOParamD["_EVreal"] = tfa::allReduce(evr, MPI_MAX);
    sim.IOParamD["_EVimag"] = tfa::allReduce(evi, MPI_MAX);
}

TF_Func void computeRealEV_GershgorinMat(tfa::Simulation &sim)
{
    // NOTE: since this depends on the geometry only, it can be calculated once
    // at the start of the simulation.

    auto &ux_N = tfa::getField(sim, "ux_N");
    uint32_t dim = ux_N.dim;

    auto &vol = tfa::meshCellVolume(sim, dim);
    auto &af  = tfa::meshFaceArea(sim, dim);
    auto &dn  = faceCellDist(sim);
    auto &fim = tfa::meshFaceInnerMask(sim, dim);

    auto &SF  = tfa::getMatrix(sim, "SumFaces");

    auto tmp = tfa::getTmpFieldS(sim, af);
    auto gg  = tfa::getTmpFieldS(sim, SF, *tmp);

    // tmp_F = (A/|d*n|)_F
    tfa::oper_copy(af, *tmp);
    tfa::oper_div(*tmp, dn);

    // If we do not want the boundary faces to contribute to the sum for each
    // cell, we just set the values of 'tmp' to zero at the boundary faces. We
    // achieve this by multiplying 'tmp' by 'fim' (face inner mask), where
    // 'fim' is 1 for inner faces (faces with exactly two neighbour cells) and
    // 0 for boundary faces.
    tfa::oper_prod(fim, *tmp);

    // gg_C = (sum(faces) tmp_f)/vol_C
    tfa::oper_prod(SF, *tmp, *gg);
    tfa::oper_div(*gg, vol);

    double kinVisc = sim.IOParamD["kinVisc"];
    sim.IOParamD["_EVreal"] = kinVisc*tfa::max(tfa::oper_max(*gg));
}

TF_Func bool computeImagEV_GershgorinMat(tfa::Simulation &sim)
{
    // NOTE: since this depends on the velocity field, it has to be updated
    // every iteration.

    auto &ux_f = tfa::getField(sim, "ux_F");
    auto &uy_f = tfa::getField(sim, "uy_F");
    auto &uz_f = tfa::getField(sim, "uz_F");

    // In a real simulation using the fractional step method, we will have the
    // velocity at the faces already, so we just get it. There is no need to
    // interpolate the velocity from the nodes.
    // auto &ux_f = tfa::getField(sim, "ux_F");
    // auto &uy_f = tfa::getField(sim, "uy_F");
    // auto &uz_f = tfa::getField(sim, "uz_F");

    uint32_t dim = ux_f.dim;

    auto &nx  = tfa::meshFaceNx(sim, dim);
    auto &ny  = tfa::meshFaceNy(sim, dim);
    auto &nz  = tfa::meshFaceNz(sim, dim);
    auto &af  = tfa::meshFaceArea(sim, dim);
    auto &fim = tfa::meshFaceInnerMask(sim, dim);
    auto &vol = tfa::meshCellVolume(sim, dim);

    auto &SF  = tfa::getMatrix(sim, "SumFaces");
    auto tmp = tfa::getTmpFieldS(sim, ux_f);
    auto gg  = tfa::getTmpFieldS(sim, SF, *tmp);

    // tmp_F = (A*|u·n|)_F
    tfa::oper_fmadd(ux_f, nx, *tmp, 1.0, 0.0);
    tfa::oper_fmadd(uy_f, ny, *tmp, 1.0, 1.0);
    tfa::oper_fmadd(uz_f, nz, *tmp, 1.0, 1.0);
    tfa::oper_max(*tmp, *tmp, 1.0, -1.0);
    tfa::oper_prod(af, *tmp);

    // If we do not want the boundary faces to contribute to the sum for each
    // cell, we just set the values of 'tmp' to zero at the boundary faces. We
    // achieve this by multiplying 'tmp' by 'fim' (face inner mask), where
    // 'fim' is 1 for inner faces (faces with exactly two neighbour cells) and
    // 0 for boundary faces.
    tfa::oper_prod(fim, *tmp);

    // gg_C = (sum(faces) tmp_f)/vol_C
    tfa::oper_prod(SF, *tmp, *gg);
    tfa::oper_div(*gg, vol);

    sim.IOParamD["_EVimag"] = 0.5*tfa::max(tfa::oper_max(*gg));
    return tfa::Iter_Continue;
}

tfa::Field &computeLambdaTilde(tfa::Simulation &sim)
{
    auto dim = tfa::getField(sim,"ux_N").dim;
    auto &result = tfa::getOrCreateField(sim, dim, "Lambda", "Faces");
    
    auto &af = tfa::meshFaceArea(sim,dim);
    auto &vol = calcFaceStaggeredVolumes(sim);
    

    tfa::oper_axpy(af,result,1.0,0.0);
    tfa::oper_prod(af,result,sim.IOParamD["kinVisc"]);
    tfa::oper_div(result,vol);
    
    return result;
}

tfa::Field &computeAbsMassFluxes(tfa::Simulation &sim)
{
    auto &ux_f = tfa::getField(sim, "ux_F");
    auto &uy_f = tfa::getField(sim, "uy_F");
    auto &uz_f = tfa::getField(sim, "uz_F");
    auto dim = ux_f.dim;

    auto &nx  = tfa::meshFaceNx(sim,dim);
    auto &ny  = tfa::meshFaceNy(sim,dim);
    auto &nz  = tfa::meshFaceNz(sim,dim);
    auto &af  = tfa::meshFaceArea(sim,dim);

    auto &result = tfa::getOrCreateField(sim, dim, "Fs", "Faces");
    
    tfa::oper_fmadd(ux_f, nx, result, 1.0, 0.0);
    tfa::oper_fmadd(uy_f, ny, result, 1.0, 1.0);
    tfa::oper_fmadd(uz_f, nz, result, 1.0, 1.0);
    tfa::oper_prod(af,result);
    tfa::oper_apply(result,[] (double x) {return fabs(x);});

    return result;    
}

double compute_AlgEigCD(tfa::Field &Fs, double factor, tfa::Simulation &sim)
{
    uint32_t resultSize = tfa::getNumEntries(Fs);
    double eigenbound = 0.0; 
      tfa::IJK F,C;
      const auto &m = tfa::getSMesh(sim);
	    auto offset = (uint32_t)m.topo.firstOwnedFaces[tfa::mpiRank()];
      double contr, c_contr;
      for (uint32_t it = 0; it < resultSize; it++)
      {
        F = tfa::getFaceIJKFromId(m,it);
        std::vector<tfa::IJK> cnb = getFaceNbCells(F,m);
        tfa::IJK CNB;
        contr=0.0;
        for(long unsigned int c=0; c<cnb.size(); c++)
        {
          CNB = cnb.at(c);
          c_contr=0.0;
          for(int f=0; f<6; ++f)
          {
            uint32_t ff = (uint32_t)tfa::getFaceIndex(m,CNB.i,CNB.j,CNB.k,f)-offset;
            tfa::IJK FF = tfa::getFaceIJKFromId(m,(long int)ff);
            if(!tfa::isBoundaryFace(m,FF))
              c_contr += Fs.array[TO_U32(ff)];
          }
          contr += c_contr/tfa::calcCellVolume(m,CNB);
        }
        eigenbound = std::max(eigenbound,contr);
      }
    return factor*eigenbound;    
}

TF_Func void computeImagEV_AlgEigCD(tfa::Simulation &sim)
{
  auto &Fs = computeAbsMassFluxes(sim);
  sim.IOParamD["_EVimag"] = compute_AlgEigCD(Fs,0.25,sim);
}

TF_Func void computeRealEV_AlgEigCD(tfa::Simulation &sim)
{
  auto &Fs = computeLambdaTilde(sim);
  sim.IOParamD["_EVreal"] = compute_AlgEigCD(Fs,1.0,sim);
}

TF_Func void computeEV_AlgEigCD(tfa::Simulation &sim)
{
  computeRealEV_AlgEigCD(sim);
  computeImagEV_AlgEigCD(sim);
}

TF_Func void printEV(tfa::Simulation &sim)
{
  tfa::info("Eigenvalues:\n");
  tfa::info("\tReal = %e\n",sim.IOParamD["_EVreal"]);
  tfa::info("\tImag = %e\n",sim.IOParamD["_EVimag"]);
}
