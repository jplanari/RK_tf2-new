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

extern "C"
{

std::vector<double> calcFacesYPos(tfa::Simulation &sim)
{
    const auto &mesh  = tfa::getSMesh(sim);
    std::set<double, fLess> facesY;
    uint32_t NfY = TO_U32(sim.IOParamI["NY"]+1);
    const auto faces = tfa::getNumFaces(mesh);

    for (int f = 0; f < faces; ++f)
    {
        const auto &F = tfa::getFaceIJKFromId(mesh,f);
	double fy = tfa::calcFaceCentroid(mesh,F)[1];
        double nvy = tfa::calcFaceNormal(mesh,F)[1];
        if (fabs(nvy) > 0.99)
        {
            facesY.insert(fy);
            if (facesY.size() == NfY)
            {
                break;
            }
        }
    }
    TF_uAssert(facesY.size() <= NfY);

    std::vector<double> localFY;
    localFY.reserve(facesY.size());
    for (auto it : facesY)
    {
        localFY.push_back(it);
    }
    auto tmp = tfa::allGatherV(localFY);
    std::sort(tmp.begin(), tmp.end());

    std::vector<double> finalFY;
    finalFY.reserve(NfY);
    auto p = tmp.begin();
    finalFY.push_back(*p);
    p++;
    for (auto pe = tmp.end(); p != pe; p++)
    {
        if (fabs(*p-finalFY.back()) > snap)
        {
            finalFY.push_back(*p);
        }
    }

    TF_uAssert(finalFY.size() == NfY, "my size is "<<finalFY.size());
    return finalFY;
}

/*void initAvgField(tfa::Simulation &sim)
{
    auto &fab = tfa::getOrCreateField(sim, sim.IOParamI["nSims"], "avgBottom", "Faces");
    const auto &mesh = tfa::getSMesh(sim);
    std::vector<double> buf_fab(tfa::getNumFaces(mesh));
    const auto bound_faces = tfa::getNumBFaces(mesh);
    for (int bf = 0; bf<bound_faces; ++bf)
    {
        const auto &BF = tfa::getBFaceIJKFromId(mesh,bf);
        if (tfa::getBoCoId(mesh,BF) == tfa::SMesh::Bottom)
        {
            buf_fab[bf] = tfa::calcFaceArea(mesh,BF);
        }
    }
    tfa::oper_setData(fab, buf_fab.data());
    tfa::Field one(sim, sim.IOParamI["nSims"], "Faces");
    tfa::oper_setConst(one, 1.0);
    sim.IOParamD["bottomArea"] = tfa::oper_dot(one, fab)[0];
}*/


void initAvgField(tfa::Simulation &sim)
{
    auto &fab = tfa::getOrCreateField(sim, 1, "avgBottom", "Faces");
    const auto &mesh = tfa::getSMesh(sim);
    std::vector<double> buf_fab(tfa::getNumEntries(fab));
    const auto faces = tfa::getNumFaces(mesh);
    for (uint32_t it = 0; it < faces; it++)
    {
        const auto &F = tfa::getFaceIJKFromId(mesh, it);
        if (tfa::getBoCoId(mesh, F) == 3)
        {
            const auto &areas = calcCellFacesArea(mesh, F);
            buf_fab[it] = areas[1];
            // areas[1] because boco 3 is an 'y' face; for 'x' faces it would
            // have been areas[0], for 'z' faces it would have been areas[2].
        }
    }
    tfa::oper_setData(fab, buf_fab.data());
    tfa::Field one(sim, 1, "Faces");
    tfa::oper_setConst(one, 1.0);
    sim.IOParamD["bottomArea"] = tfa::oper_dot(one, fab)[0];
}

void load(tfa::Simulation &sim)
{
    const auto &name = tfa::substituteParams(sim, sim.IOParamS["load"]);
    tfa::loadFields(sim, name);
}

void averageXZ(tfa::Simulation &sim)
{
    const auto &facesYPos = calcFacesYPos(sim);
    auto NY = facesYPos.size()-1;
    std::vector<double> avguL(NY);
    std::vector<double> avgwL(NY);
    std::vector<double> avgpL(NY);
    std::vector<double> volumesL(NY);
    std::vector<double> yL(NY, -DBL_MAX);
    std::vector<double> avgdudyL(NY+1);
    std::vector<double> avgdwdyL(NY+1);
    std::vector<double> areasL(NY+1);
    std::vector<double> yfL(NY+1, -DBL_MAX);

    auto &GY = tfa::getMatrix(sim, "GradY_NF");
    auto &ux = tfa::getField(sim, "ens_avg(ux)");
    auto &uz = tfa::getField(sim, "ens_avg(uz)");
    auto &dudy = tfa::getOrCreateField(sim, 1, "dudy", "Faces");
    tfa::Field dwdy(sim, 1, "Faces");

    tfa::oper_prod(GY, ux, dudy);
    tfa::oper_prod(GY, uz, dwdy);
    auto &fab = tfa::getField(sim, "avgBottom");
    double bottomArea = sim.IOParamD["bottomArea"];
    double sdudy = fabs(tfa::oper_dot(fab, dudy)[0])/bottomArea;
    tfa::info("Calculated Re_tau is %.5e\n",sdudy);

    const auto &mesh = tfa::getSMesh(sim);
    std::vector<double> bufferu(tfa::getNumNodes(mesh));
    std::vector<double> bufferw(tfa::getNumNodes(mesh));
    std::vector<double> bufferp(tfa::getNumNodes(mesh));
    tfa::oper_getData(ux, bufferu.data());
    tfa::oper_getData(uz, bufferw.data());
    auto cells = tfa::getNumCells(mesh);
    for (uint32_t c = 0; c < cells; c++)
    {
        const auto &C = tfa::getCellIJKFromId(mesh,c);
	double cy = tfa::calcCellCentroid(mesh,C)[1];
        double vol = tfa::calcCellVolume(mesh,C);
        uint64_t idp = bsearchVectorE(facesYPos, cy, snap);
        avguL[idp] += vol*bufferu[c];
        avgwL[idp] += vol*bufferw[c];
        volumesL[idp] += vol;
        yL[idp] = cy;
    }

    const auto &outP = sim.IOParamS["outPrefix"];
    const auto &avgu = allReduce(avguL, MPI_SUM);
    const auto &avgw = allReduce(avgwL, MPI_SUM);
    const auto &volumes = allReduce(volumesL, MPI_SUM);
    const auto &y = allReduce(yL, MPI_MAX);
    const auto Re = sim.IOParamD["Re_tau"];
    const auto nsims = sim.IOParamI["nSims"];
    const auto fct = sim.IOParamD["RKfct"];
    const auto method = sim.IOParamS["RKmethod"];
    {
    std::ofstream out("results/profile_Re"+std::to_string(Re)+"_fct"+std::to_string(fct)+"_nSims"+std::to_string(nsims)+"_cells_"+method+".dat");
    TF_uAssert(out.is_open(), "could not open filename for output");
    out<<"#y y+ uavg wavg pavg\n";
    out<<std::scientific<<std::setprecision(8);
    if (tfa::mpiRank() == 0)
    {
        for (uint32_t k = 0; k < NY; k++)
        {
            out<<y[k]
               <<" "<<y[k]*sdudy
               <<" "<<avgu[k]/volumes[k]
               <<" "<<avgw[k]/volumes[k]
               <<"\n";
        }
    }
    out.close();
    }

    std::vector<double> bufferdudy(tfa::getNumFaces(mesh));
    std::vector<double> bufferdwdy(tfa::getNumFaces(mesh));
    tfa::oper_getData(dudy, bufferdudy.data());
    tfa::oper_getData(dwdy, bufferdwdy.data());
    auto faces = tfa::getNumFaces(mesh);
    for (uint32_t f=0; f<faces; ++f)
    {
        const auto &F = tfa::getFaceIJKFromId(mesh,f);
        auto nvy = tfa::calcFaceNormal(mesh,F)[1];
        if (fabs(nvy) > 0.99)
        {
            double cy = tfa::calcFaceCentroid(mesh,F)[1];
	    const auto &areas = calcCellFacesArea(mesh,F);
            double area = fabs(areas*tfa::calcFaceNormal(mesh,F));
            uint64_t idp = bsearchVectorE(facesYPos, cy, snap);
            avgdudyL[idp] += area*bufferdudy[f];
            avgdwdyL[idp] += area*bufferdwdy[f];
            areasL[idp] += area;
            yfL[idp] = cy;
        }
    }

    const auto &avgdudy = allReduce(avgdudyL, MPI_SUM);
    const auto &avgdwdy = allReduce(avgdwdyL, MPI_SUM);
    const auto &areas = allReduce(areasL, MPI_SUM);
    const auto &yf = allReduce(yfL, MPI_MAX);
    {
    std::ofstream out("results/profile_Re"+std::to_string(Re)+"_fct"+std::to_string(fct)+"_nSims"+std::to_string(nsims)+"_faces_"+method+".dat");
    TF_uAssert(out.is_open(), "could not open filename for output");
    out<<"#y y+ duavg/dy dwavg/dy\n";
    out<<std::scientific<<std::setprecision(8);
    if (tfa::mpiRank() == 0)
    {
        for (uint32_t k = 0; k < NY+1; k++)
        {
            out<<yf[k]
               <<" "<<yf[k]*sdudy
               <<" "<<avgdudy[k]/areas[k]
               <<" "<<avgdwdy[k]/areas[k]
               <<"\n";
        }
    }
    out.close();
    }
}

void ReynoldsStresses(tfa::Simulation &sim)
{
    auto &avg_u = tfa::getField(sim, "ens_avg(ux)");
    auto &avg_v = tfa::getField(sim, "ens_avg(uy)");
    auto &avg_w = tfa::getField(sim, "ens_avg(uz)");
    auto &avg_uu = tfa::getField(sim, "ens_avg(ux*ux)");
    auto &avg_uv = tfa::getField(sim, "ens_avg(ux*uy)");
    auto &avg_uw = tfa::getField(sim, "ens_avg(ux*uz)");
    auto &avg_vv = tfa::getField(sim, "ens_avg(uy*uy)");
    auto &avg_vw = tfa::getField(sim, "ens_avg(uy*uz)");
    auto &avg_ww = tfa::getField(sim, "ens_avg(uz*uz)");

    tfa::Field R_uu(sim, 1, "Nodes");
    tfa::Field R_uv(sim, 1, "Nodes");
    tfa::Field R_uw(sim, 1, "Nodes");
    tfa::Field R_vv(sim, 1, "Nodes");
    tfa::Field R_vw(sim, 1, "Nodes");
    tfa::Field R_ww(sim, 1, "Nodes");

    tfa::oper_fmadd(avg_u, avg_u, R_uu, 1.0, 0.0);
    tfa::oper_axpy(avg_uu, R_uu, 1.0, -1.0);

    tfa::oper_fmadd(avg_u, avg_v, R_uv, 1.0, 0.0);
    tfa::oper_axpy(avg_uv, R_uv, 1.0, -1.0);

    tfa::oper_fmadd(avg_u, avg_w, R_uw, 1.0, 0.0);
    tfa::oper_axpy(avg_uw, R_uw, 1.0, -1.0);

    tfa::oper_fmadd(avg_v, avg_v, R_vv, 1.0, 0.0);
    tfa::oper_axpy(avg_vv, R_vv, 1.0, -1.0);

    tfa::oper_fmadd(avg_v, avg_w, R_vw, 1.0, 0.0);
    tfa::oper_axpy(avg_vw, R_vw, 1.0, -1.0);

    tfa::oper_fmadd(avg_w, avg_w, R_ww, 1.0, 0.0);
    tfa::oper_axpy(avg_ww, R_ww, 1.0, -1.0);

    const auto &facesYPos = calcFacesYPos(sim);
    auto NY = facesYPos.size()-1;
    std::vector<double> avgLR_uu(NY);
    std::vector<double> avgLR_uv(NY);
    std::vector<double> avgLR_uw(NY);
    std::vector<double> avgLR_vv(NY);
    std::vector<double> avgLR_vw(NY);
    std::vector<double> avgLR_ww(NY);
    std::vector<double> volumesL(NY);
    std::vector<double> yL(NY, -DBL_MAX);
    const auto &mesh = tfa::getSMesh(sim);
    std::vector<double> buffer_uu(tfa::getNumNodes(mesh));
    std::vector<double> buffer_uv(tfa::getNumNodes(mesh));
    std::vector<double> buffer_uw(tfa::getNumNodes(mesh));
    std::vector<double> buffer_vv(tfa::getNumNodes(mesh));
    std::vector<double> buffer_vw(tfa::getNumNodes(mesh));
    std::vector<double> buffer_ww(tfa::getNumNodes(mesh));
    tfa::oper_getData(R_uu, buffer_uu.data());
    tfa::oper_getData(R_uv, buffer_uv.data());
    tfa::oper_getData(R_uw, buffer_uw.data());
    tfa::oper_getData(R_vv, buffer_vv.data());
    tfa::oper_getData(R_vw, buffer_vw.data());
    tfa::oper_getData(R_ww, buffer_ww.data());
    auto const cells = tfa::getNumCells(mesh);
    for (int c = 0; c<cells; ++c)
    {
        const auto &C = tfa::getCellIJKFromId(mesh,c);
        double cy = tfa::calcCellCentroid(mesh,C)[1];
        double vol = tfa::calcCellVolume(mesh,C);
        uint64_t idp = bsearchVectorE(facesYPos, cy, snap);
        avgLR_uu[idp] += vol*buffer_uu[c];
        avgLR_uv[idp] += vol*buffer_uv[c];
        avgLR_uw[idp] += vol*buffer_uw[c];
        avgLR_vv[idp] += vol*buffer_vv[c];
        avgLR_vw[idp] += vol*buffer_vw[c];
        avgLR_ww[idp] += vol*buffer_ww[c];
        volumesL[idp] += vol;
        yL[idp] = cy;
    }

    auto &GY = tfa::getMatrix(sim, "GradY_NF");
    tfa::Field dudy(sim, 1, "Faces");
    tfa::oper_prod(GY, avg_u, dudy);
    auto &fab = tfa::getField(sim, "avgBottom");
    double bottomArea = sim.IOParamD["bottomArea"];
    double sdudy = fabs(tfa::oper_dot(fab, dudy)[0])/bottomArea;

    const auto &outP = sim.IOParamS["outPrefix"];
    const auto &avgR_uu = allReduce(avgLR_uu, MPI_SUM);
    const auto &avgR_uv = allReduce(avgLR_uv, MPI_SUM);
    const auto &avgR_uw = allReduce(avgLR_uw, MPI_SUM);
    const auto &avgR_vv = allReduce(avgLR_vv, MPI_SUM);
    const auto &avgR_vw = allReduce(avgLR_vw, MPI_SUM);
    const auto &avgR_ww = allReduce(avgLR_ww, MPI_SUM);
    const auto &volumes = allReduce(volumesL, MPI_SUM);
    const auto Re = sim.IOParamD["Re_tau"];
    const auto fct = sim.IOParamD["RKfct"];
    const auto nsims = sim.IOParamI["nSims"];
    const auto method = sim.IOParamS["RKmethod"];
    const auto &y = allReduce(yL, MPI_MAX);
    {
    std::ofstream out("results/reyStress_Re"+std::to_string(Re)+"_fct"+std::to_string(fct)+"_nSims_"+std::to_string(nsims)+"_faces_"+method+".dat");
    TF_uAssert(out.is_open(), "could not open filename for output");
    out<<"#y y+ R_uu R_vv R_ww R_uv R_uw R_vw\n";
    out<<std::scientific<<std::setprecision(8);
    if (tfa::mpiRank() == 0)
    {
        for (uint32_t k = 0; k < NY; k++)
        {
            out<<y[k]
               <<" "<<y[k]*sdudy
               <<" "<<avgR_uu[k]/volumes[k]
               <<" "<<avgR_vv[k]/volumes[k]
               <<" "<<avgR_ww[k]/volumes[k]
               <<" "<<avgR_uv[k]/volumes[k]
               <<" "<<avgR_uw[k]/volumes[k]
               <<" "<<avgR_vw[k]/volumes[k]
               <<"\n";
        }
    }
    out.close();
    }
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
    auto &uxux = tfa::getField(sim, "avg(ux_N*ux_N)");
    auto &uxuy = tfa::getField(sim, "avg(ux_N*uy_N)");
    auto &uxuz = tfa::getField(sim, "avg(ux_N*uz_N)");
    auto &uyuy = tfa::getField(sim, "avg(uy_N*uy_N)");
    auto &uyuz = tfa::getField(sim, "avg(uy_N*uz_N)");
    auto &uzuz = tfa::getField(sim, "avg(uz_N*uz_N)");

    auto &ensAvg_ux = tfa::getOrCreateField(sim, 1, "ens_avg(ux)", "Nodes");
    auto &ensAvg_uy = tfa::getOrCreateField(sim, 1, "ens_avg(uy)", "Nodes");
    auto &ensAvg_uz = tfa::getOrCreateField(sim, 1, "ens_avg(uz)", "Nodes");
    auto &ensAvg_uxux = tfa::getOrCreateField(sim, 1, "ens_avg(ux*ux)", "Nodes");
    auto &ensAvg_uxuy = tfa::getOrCreateField(sim, 1, "ens_avg(ux*uy)", "Nodes");
    auto &ensAvg_uxuz = tfa::getOrCreateField(sim, 1, "ens_avg(ux*uz)", "Nodes");
    auto &ensAvg_uyuy = tfa::getOrCreateField(sim, 1, "ens_avg(uy*uy)", "Nodes");
    auto &ensAvg_uyuz = tfa::getOrCreateField(sim, 1, "ens_avg(uy*uz)", "Nodes");
    auto &ensAvg_uzuz = tfa::getOrCreateField(sim, 1, "ens_avg(uz*uz)", "Nodes");

    apply_ensemble_average(ux, ensAvg_ux);
    apply_ensemble_average(uy, ensAvg_uy);
    apply_ensemble_average(uz, ensAvg_uz);
    apply_ensemble_average(uxux, ensAvg_uxux);
    apply_ensemble_average(uxuy, ensAvg_uxuy);
    apply_ensemble_average(uxuz, ensAvg_uxuz);
    apply_ensemble_average(uyuy, ensAvg_uyuy);
    apply_ensemble_average(uyuz, ensAvg_uyuz);
    apply_ensemble_average(uzuz, ensAvg_uzuz);
}

}
