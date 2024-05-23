#include <algorithm>
#include <cfloat>
#include <iomanip>
#include "tfa/DataIO.h"
#include "tfa/IO.h"
#include "tfa/Simulation.h"
#include "tfa/Opers.h"
#include "tfa/ll_Opers.h"

static const double snap = 1e-8;

struct fLess
{
    constexpr bool operator()(const double &lhs, const double &rhs) const
    {
        return (lhs < rhs && fabs(lhs-rhs) > snap);
    }
};

uint64_t bsearchVectorE(const std::vector<double> &vec, double target, double eps)
{
    // Given a vector 'vec' such that vec[i]>vec[i-1], find the largest i such
    // that vec[i]<=target using a binary search.
    uint64_t li = 0, ui = vec.size()-1;
    while ((ui-li)>1)
    {
        uint64_t mi = (li+ui)/2;
        if (fabs(vec[mi]-target) < eps)
        {
            return mi;
        }
        else if (vec[mi]>target)
        {
            ui = mi;
        }
        else
        {
            li = mi;
        }
    }
    return vec[ui] <= target? ui : li;
}

std::vector<double> allReduce(const std::vector<double> &in, MPI_Op op, MPI_Comm w = MPI_COMM_WORLD)
{
    std::vector<double> result(in.size());
    int32_t sz = TO_S32(in.size());
    tfa::checkr(MPI_Allreduce(in.data(), result.data(), sz, MPI_DOUBLE, op, w),
                "allreduce");
    return result;
}

TF_Func void SetUp_Matrices(tfa::Simulation &sim)
{
    // This function exists basically to initialize all the fields required by
    // the official FSM, with the correct number of components, i.e., one
    // component per simulation.
    const auto &meshName = tfa::getMeshName(sim);

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
}

extern "C"
{

std::vector<double> calcFacesZPos(tfa::Simulation &sim)
{
    const auto &mesh  = tfa::getSMesh(sim);
    std::set<double, fLess> facesZ;
    uint32_t NfZ = TO_U32(sim.IOParamI["NZ"]+1);
    const auto faces = tfa::getNumFaces(mesh);

    for (int f = 0; f < faces; ++f)
    {
        const auto &F = tfa::getFaceIJKFromId(mesh,f);
	      double fz = tfa::calcFaceCentroid(mesh,F)[2];
        double nvz = tfa::calcFaceNormal(mesh,F)[2];
        if (fabs(nvz) > 0.99)
        {
            facesZ.insert(fz);
            if (facesZ.size() == NfZ)
            {
                break;
            }
        }
    }
    TF_uAssert(facesZ.size() <= NfZ);

    std::vector<double> localFZ;
    localFZ.reserve(facesZ.size());
    for (auto it : facesZ)
    {
        localFZ.push_back(it);
    }
    auto tmp = tfa::allGatherV(localFZ);
    std::sort(tmp.begin(), tmp.end());

    std::vector<double> finalFZ;
    finalFZ.reserve(NfZ);
    auto p = tmp.begin();
    finalFZ.push_back(*p);
    p++;
    for (auto pe = tmp.end(); p != pe; p++)
    {
        if (fabs(*p-finalFZ.back()) > snap)
        {
            finalFZ.push_back(*p);
        }
    }

    TF_uAssert(finalFZ.size() == NfZ, "my size is "<<finalFZ.size());
    return finalFZ;
}

void initAvgField(tfa::Simulation &sim)
{
    auto &fab = tfa::getOrCreateField(sim, 1, "avgBack", "Faces");
    const auto &mesh = tfa::getSMesh(sim);
    std::vector<double> buf_fab(tfa::getNumEntries(fab));
    const auto faces = tfa::getNumFaces(mesh);
    for (uint32_t it = 0; it < faces; it++)
    {
        const auto &F = tfa::getFaceIJKFromId(mesh, it);
        if (tfa::getBoCoId(mesh, F) == 4)
        {
            const auto &areas = calcCellFacesArea(mesh, F);
            buf_fab[it] = areas[2];
            // areas[1] because boco 3 is an 'y' face; for 'x' faces it would
            // have been areas[0], for 'z' faces it would have been areas[2].
        }
    }
    tfa::oper_setData(fab, buf_fab.data());
    tfa::Field one(sim, 1, "Faces");
    tfa::oper_setConst(one, 1.0);
    sim.IOParamD["backArea"] = tfa::oper_dot(one, fab)[0];
}



void load(tfa::Simulation &sim)
{
    const auto &name = tfa::substituteParams(sim, sim.IOParamS["load"]);
    tfa::loadFields(sim, name);
}

void apply_ensemble_average(tfa::Field &src, tfa::Field &dest)
{
    TF_uAssert(dest.dim == 1);
    TF_uAssert(&tfa::getDomain(src) == &tfa::getDomain(dest));

    std::vector<double> srcBuffer(tfa::getNumEntries(src));
    std::vector<double> destBuffer(tfa::getNumEntries(dest));

    tfa::oper_getData(src, srcBuffer.data());

    auto srcElem = tfa::getDomain(src).localSize;
    for (uint32_t k = 0; k < srcElem; k++)
    {
        double avg = 0.0;
        for (uint32_t d = 0; d < src.dim; d++) avg += srcBuffer[k*src.dim+d];
        destBuffer[k] = avg/TO_F64(src.dim);
    }

    tfa::oper_setData(dest, destBuffer.data());
}

void ensemble_avg(tfa::Simulation &sim)
{
    auto &ux = tfa::getField(sim, "avg(ux_N)");
    auto &uy = tfa::getField(sim, "avg(uy_N)");
    auto &uz = tfa::getField(sim, "avg(uz_N)");
    auto &T = tfa::getField(sim, "avg(T_N)");

    auto &ensAvg_ux = tfa::getOrCreateField(sim, 1, "ens_avg(ux)", "Nodes");
    auto &ensAvg_uy = tfa::getOrCreateField(sim, 1, "ens_avg(uy)", "Nodes");
    auto &ensAvg_uz = tfa::getOrCreateField(sim, 1, "ens_avg(uz)", "Nodes");
    auto &ensAvg_T = tfa::getOrCreateField(sim, 1, "ens_avg(T)", "Nodes");

    apply_ensemble_average(ux, ensAvg_ux);
    apply_ensemble_average(uy, ensAvg_uy);
    apply_ensemble_average(uz, ensAvg_uz);
    apply_ensemble_average(T, ensAvg_T);
}

void computeNusselt(tfa::Simulation &sim)
{
  auto &T = tfa::getField(sim, "ens_avg(T)");
  auto &GX = tfa::getMatrix(sim, "GradX_CC");
  auto &ID_NC = tfa::getMatrix(sim, "Id_NC");

  auto NX = sim.IOParamI["NX"];
  auto Tc = tfa::getTmpFieldS(sim, 1, "Cells");
  tfa::oper_prod(ID_NC, T, *Tc);

  auto dTdx = tfa::getTmpFieldS(sim, *Tc);
  tfa::oper_prod(GX, *Tc, *dTdx);
  
  const auto &m = tfa::getSMesh(sim);
  std::vector<double> dTdx_buf(tfa::getNumEntries(*dTdx));
  std::vector<double> Nux_buf(tfa::getNumEntries(*dTdx));
  std::vector<double> left_buf(tfa::getNumEntries(*dTdx));
  std::vector<double> right_buf(tfa::getNumEntries(*dTdx));
  tfa::oper_getData(*dTdx, dTdx_buf.data());

  auto srcElem = tfa::getDomain(*dTdx).localSize;
  double dy,dz;
  for (uint32_t k=0; k < srcElem; k++)
  {
    auto C = tfa::getCellIJKFromId(m,k);
    tfa::info("cell %d, C.i=%d, NX=%d\n",k,C.i,NX);
    auto Ftop = tfa::getFaceIJKFromId(m,tfa::getFaceIndex(m,C.i,C.j,C.k,3));
    auto Fbot = tfa::getFaceIJKFromId(m,tfa::getFaceIndex(m,C.i,C.j,C.k,2));

    dy = tfa::calcFaceCentroid(m,Ftop)[1]-tfa::calcFaceCentroid(m,Fbot)[1];

    auto Ffro = tfa::getFaceIJKFromId(m,tfa::getFaceIndex(m,C.i,C.j,C.k,5));
    auto Fbac = tfa::getFaceIJKFromId(m,tfa::getFaceIndex(m,C.i,C.j,C.k,4));

    dz = tfa::calcFaceCentroid(m,Ffro)[2] - tfa::calcFaceCentroid(m,Fbac)[2];

    Nux_buf[k] = dTdx_buf[k]*dy*dz;

    if (C.i==0) left_buf[k] = -1.0;
    else left_buf[k] = 0.0;
    
    if (C.i==NX-1) right_buf[k] = 1.0;
    else right_buf[k] = 0.0;
  }
  
  auto Nu = tfa::getTmpFieldS(sim, *dTdx);
  tfa::oper_setData(*Nu, Nux_buf.data());
  
  auto left = tfa::getTmpFieldS(sim, *Nu);
  tfa::oper_setData(*left,left_buf.data());
  tfa::info("Nusselt = %e\n",tfa::oper_dot(*Nu,*left)[0]/4);
}
}
