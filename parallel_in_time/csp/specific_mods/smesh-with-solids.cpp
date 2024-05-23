#include "tfa/Opers.h"
#include "tfa/Simulation.h"

const auto &CI = tfa::getCellIndex;

static double cavWidth  = 0.025;
static double cavHeight = 0.050;
static double cavDepth  = 0.025;

static double calcBndDist(const tfa::SMesh &m, tfa::IJK N)
{
    auto nvc = tfa::getSize(m);
    double result = 0.0;
    if      (N.i == -1)     result = tfa::calcNodeDX(m, 0, -1);
    else if (N.i == nvc[0]) result = tfa::calcNodeDX(m, nvc[0], nvc[0]-1);
    else if (N.j == -1)     result = tfa::calcNodeDY(m, 0, -1);
    else if (N.j == nvc[1]) result = tfa::calcNodeDY(m, nvc[1], nvc[1]-1);
    else if (N.k == -1)     result = tfa::calcNodeDZ(m, 0, -1);
    else if (N.k == nvc[2]) result = tfa::calcNodeDZ(m, nvc[2], nvc[2]-1);
    else    TF_uAssert(false, "BUG! not a boundary node");

    return result;
}

bool isSolid(const tfa::SMesh &m, int i, int j, int k)
{
    tfa::IJK C{i, j, k};
    tfa::V3 cen = tfa::calcCellCentroid(m, C);

    // no solids
    // return false;

    // return (cen[0] <= 0.5 && cen[1] <= 0.7);
    // return (cen[0] >= 0.25 && cen[0] <= 0.75 && cen[1] >= 0.25 && cen[1] <= 0.75);

    // return (cen[0] >= 0.25 && cen[0] <= 0.75 && cen[1] >= 1.20 && cen[1] <= 1.60);

    // return (cen[0] <= 0.3 && cen[1] >= 1.20 && cen[1] <= 1.30);
    // return (cen[0] <= 0.15 && cen[1] >= 1.20);

    // 4 solids
    // return ((cen[0] >= 0.45 && cen[0] <= 0.55 && cen[1] >= 0.45 && cen[1] <= 0.55) ||
    //         (cen[0] >= 0.95 && cen[0] <= 1.05 && cen[1] >= 0.60 && cen[1] <= 0.70) ||
    //         (cen[0] >= 0.95 && cen[0] <= 1.05 && cen[1] >= 0.30 && cen[1] <= 0.40) ||
    //         (cen[0] >= 1.45 && cen[0] <= 1.55 && cen[1] >= 0.45 && cen[1] <= 0.55));

    // csp-2d
    // double domHeight = tfa::getLength(m)[1];
    // return cen[0] < cavWidth &&
    //        (cen[1] < 0.5*(domHeight-cavHeight) || cen[1] > 0.5*(domHeight+cavHeight));

    // csp-3d
    double domHeight = tfa::getLength(m)[1];
    double domDepth  = tfa::getLength(m)[2];
    bool inFluid = (cen[0] >= cavWidth ||
                   (cen[1] >= 0.5*(domHeight-cavHeight) && cen[1] <= 0.5*(domHeight+cavHeight) &&
                    cen[2] >= 0.5*(domDepth-cavDepth) && cen[2] <= 0.5*(domDepth+cavDepth)));
    return not inFluid;
}

TF_Func void cav_2parts(const tfa::SMeshVertParams &p)
{
    TF_uAssert(p.pars.size() == 3, "cavXVerts, want 3 param here: "
               "cavity_width, fraction_of_nodes_in_cavity, concentration_factor");
    double cavWidth;
    double pct1;
    double alpha;
    TF_uAssert(tfa::isNumber(p.pars[0], cavWidth), "expected number, got '"+p.pars[0]+"'");
    TF_uAssert(tfa::isNumber(p.pars[1], pct1),     "expected number, got '"+p.pars[1]+"'");
    TF_uAssert(tfa::isNumber(p.pars[2], alpha),    "expected number, got '"+p.pars[2]+"'");
    int32_t nc = p.num-1;
    int32_t n1 = TO_S32(nc*pct1);
    int32_t n2 = nc-n1;
    double len1 = cavWidth;
    double len2 = p.last-p.first-cavWidth;
    for (int32_t it = 0; it < n1; it++)
    {
        p.verts[it] = p.first+it*len1/n1;
    }
    for (int32_t it = n1; it <= nc; it++)
    {
        // p.verts[it] = p.first+len1+(it-n1)*len2/n2;
        // p.verts[it] = p.first+len1+len2*(tanh(alpha*(it-n1)/n2)/tanh(alpha));
        p.verts[it] = p.first+len1+len2*(1.0-tanh(alpha*(n2-it+n1)/n2)/tanh(alpha));
    }
}

TF_Func void cav_3parts(const tfa::SMeshVertParams &p)
{
    TF_uAssert(p.pars.size() == 3, "cav_3parts, want 3 params here:"
               "cavity_dimension, fraction_of_nodes_in_cavity, "
               "concentration_factor_outside");
    double cavDim;
    double cavPct;
    double alpha;
    TF_uAssert(tfa::isNumber(p.pars[0], cavDim), "expected number, got '"+p.pars[0]+"'");
    TF_uAssert(tfa::isNumber(p.pars[1], cavPct), "expected number, got '"+p.pars[1]+"'");
    TF_uAssert(tfa::isNumber(p.pars[2], alpha),  "expected number, got '"+p.pars[2]+"'");
    int32_t nc = p.num-1;
    int32_t n2 = TO_S32(nc*cavPct);
    int32_t n1 = (nc-n2)/2;
    int32_t n3 = n1;
    double len1 = (p.last-p.first-cavDim)*0.5;
    double len2 = cavDim;
    double len3 = len1;
    for (int32_t it = 0; it < n1; it++)
    {
        p.verts[it] = p.first+len1*(tanh(alpha*it/n1)/tanh(alpha));
    }
    for (int32_t it = n1; it < n1+n2; it++)
    {
        p.verts[it] = p.first+len1+(it-n1)*len2/n2;
    }
    for (int32_t it = n1+n2; it <= nc; it++)
    {
        p.verts[it] = p.first+len1+len2+len3*(1.0-tanh(alpha*(n3-it+n1+n2)/n3)/tanh(alpha));
    }
}

TF_Func void init_geom(tfa::Simulation &sim)
{
    cavWidth  = tfa::getParamD(sim, "cavWidth");
    cavHeight = tfa::getParamD(sim, "cavHeight");
    cavDepth  = tfa::getParamD(sim, "cavDepth");
}

TF_Func void setup_marks(tfa::Simulation &sim)
{
    auto &solids = tfa::getOrCreateField(sim, 1, "Solids", "Nodes");
    std::vector<double> buffer(tfa::getNumEntries(solids), 0.0);

    const auto &m = tfa::getSMesh(sim);
    for (uint32_t it = 0; it < tfa::getNumNodes(m); it++)
    {
        auto N = getNodeIJKFromId(m, it);
        if (isSolid(m, N.i, N.j, N.k)) buffer[it] = 1.0;
    }

    tfa::oper_setData(solids, buffer.data());

    auto &fm = tfa::getOrCreateField(sim, "FluidMask", solids);
    tfa::oper_copy(solids, fm);
    tfa::oper_add(fm, -1.0);
    tfa::oper_scale(fm, -1.0);
}

TF_Func void zero_vel_at_solids(tfa::Simulation &sim)
{
    auto &fm = tfa::getField(sim, "FluidMask");
    auto &ux = tfa::getField(sim, "ux");
    auto &uy = tfa::getField(sim, "uy");
    auto &uz = tfa::getField(sim, "uz");
    tfa::oper_prod(fm, ux);
    tfa::oper_prod(fm, uy);
    tfa::oper_prod(fm, uz);
}

TF_Func uint32_t k_PEq_CC_FirstNb_SPD_w_solids(const tfa::KernelParams &kp)
{
    NUM_COMPONENTS(kp, 1);
    const auto &mesh = tfa::getSMesh(kp.sim);
    const auto &ijk = tfa::getCellIJKFromId(mesh, kp.row);
    auto i = ijk.i;
    auto j = ijk.j;
    auto k = ijk.k;
    double iVol = 1.0;
    tfa::V3 areas = tfa::calcCellFacesArea(mesh, ijk);
    double xcoef = areas[0]*iVol;
    double ycoef = areas[1]*iVol;
    double zcoef = areas[2]*iVol;
    uint32_t idx = kp.idx;
    double coef = 0.0;

    if (isSolid(mesh, i, j, k))
    {
        if (CI(mesh, i-1, j, k) != -1) kp.values[idx++] = 0.0;
        if (CI(mesh, i+1, j, k) != -1) kp.values[idx++] = 0.0;
        if (CI(mesh, i, j-1, k) != -1) kp.values[idx++] = 0.0;
        if (CI(mesh, i, j+1, k) != -1) kp.values[idx++] = 0.0;
        if (CI(mesh, i, j, k-1) != -1) kp.values[idx++] = 0.0;
        if (CI(mesh, i, j, k+1) != -1) kp.values[idx++] = 0.0;
        kp.values[idx] = -1.0;
    }
    else
    {
        if (CI(mesh, i-1, j, k) != -1) coef += (kp.values[idx++] = (not isSolid(mesh, i-1, j, k))*xcoef/tfa::calcCellDX(mesh, i,   i-1));
        if (CI(mesh, i+1, j, k) != -1) coef += (kp.values[idx++] = (not isSolid(mesh, i+1, j, k))*xcoef/tfa::calcCellDX(mesh, i+1, i));
        if (CI(mesh, i, j-1, k) != -1) coef += (kp.values[idx++] = (not isSolid(mesh, i, j-1, k))*ycoef/tfa::calcCellDY(mesh, j,   j-1));
        if (CI(mesh, i, j+1, k) != -1) coef += (kp.values[idx++] = (not isSolid(mesh, i, j+1, k))*ycoef/tfa::calcCellDY(mesh, j+1, j));
        if (CI(mesh, i, j, k-1) != -1) coef += (kp.values[idx++] = (not isSolid(mesh, i, j, k-1))*zcoef/tfa::calcCellDZ(mesh, k,   k-1));
        if (CI(mesh, i, j, k+1) != -1) coef += (kp.values[idx++] = (not isSolid(mesh, i, j, k+1))*zcoef/tfa::calcCellDZ(mesh, k+1, k));
        kp.values[idx] = -coef;
    }

    return 0;
}

TF_Func uint32_t k_PEq_NN_FirstNb_SPD_w_solids(const tfa::KernelParams &kp)
{
    NUM_COMPONENTS(kp, 1);
    const auto &mesh = tfa::getSMesh(kp.sim);
    const auto &N = tfa::getNodeIJKFromId(mesh, kp.row);
    auto i = N.i;
    auto j = N.j;
    auto k = N.k;

    if (tfa::isBoundaryFace(mesh, N))
    {
        auto C = tfa::getFaceNbCellsIJK(mesh, N).first;
        auto coef = (not isSolid(mesh, C.i, C.j, C.k)) ? tfa::calcFaceArea(mesh, N)/calcBndDist(mesh, N) : 0.0;
        kp.values[kp.idx+0] = -coef;
        kp.values[kp.idx+1] = +coef;
        return 0;
    }

    const tfa::V3 areas = tfa::calcCellFacesArea(mesh, N);
    const double xcoef = -areas[0];
    const double ycoef = -areas[1];
    const double zcoef = -areas[2];
    const auto &t = mesh.topo;

    if (isSolid(mesh, i, j, k))
    {
        kp.values[kp.idx+0] = 0.0;
        kp.values[kp.idx+1] = 0.0;
        kp.values[kp.idx+2] = 0.0;
        kp.values[kp.idx+3] = 0.0;
        kp.values[kp.idx+4] = 0.0;
        kp.values[kp.idx+5] = 0.0;
        kp.values[kp.idx+6] = 0.0;
    }
    else
    {
        // double xm = not (i-1 == -1   || isSolid(mesh, i-1, j, k))? xcoef/tfa::calcNodeDX(mesh, i, i-1) : 0.0*xcoef*xcoef*2.0;
        // double xp = not (i+1 == t.nx || isSolid(mesh, i+1, j, k))? xcoef/tfa::calcNodeDX(mesh, i+1, i) : 0.0*xcoef*xcoef*2.0;
        // double ym = not (j-1 == -1   || isSolid(mesh, i, j-1, k))? ycoef/tfa::calcNodeDY(mesh, j, j-1) : 0.0*ycoef*ycoef*2.0;
        // double yp = not (j+1 == t.ny || isSolid(mesh, i, j+1, k))? ycoef/tfa::calcNodeDY(mesh, j+1, j) : 0.0*ycoef*ycoef*2.0;
        // double zm = not (k-1 == -1   || isSolid(mesh, i, j, k-1))? zcoef/tfa::calcNodeDZ(mesh, k, k-1) : 0.0*zcoef*zcoef*2.0;
        // double zp = not (k+1 == t.nz || isSolid(mesh, i, j, k+1))? zcoef/tfa::calcNodeDZ(mesh, k+1, k) : 0.0*zcoef*zcoef*2.0;
        double xm = (i-1 == -1   || not isSolid(mesh, i-1, j, k))? xcoef/tfa::calcNodeDX(mesh, i, i-1) : 0.0*xcoef*xcoef*2.0;
        double xp = (i+1 == t.nx || not isSolid(mesh, i+1, j, k))? xcoef/tfa::calcNodeDX(mesh, i+1, i) : 0.0*xcoef*xcoef*2.0;
        double ym = (j-1 == -1   || not isSolid(mesh, i, j-1, k))? ycoef/tfa::calcNodeDY(mesh, j, j-1) : 0.0*ycoef*ycoef*2.0;
        double yp = (j+1 == t.ny || not isSolid(mesh, i, j+1, k))? ycoef/tfa::calcNodeDY(mesh, j+1, j) : 0.0*ycoef*ycoef*2.0;
        double zm = (k-1 == -1   || not isSolid(mesh, i, j, k-1))? zcoef/tfa::calcNodeDZ(mesh, k, k-1) : 0.0*zcoef*zcoef*2.0;
        double zp = (k+1 == t.nz || not isSolid(mesh, i, j, k+1))? zcoef/tfa::calcNodeDZ(mesh, k+1, k) : 0.0*zcoef*zcoef*2.0;
        kp.values[kp.idx+0] = (i-1 == -1   || not isSolid(mesh, i-1, j, k))*xm;
        kp.values[kp.idx+1] = (i+1 == t.nx || not isSolid(mesh, i+1, j, k))*xp;
        kp.values[kp.idx+2] = (j-1 == -1   || not isSolid(mesh, i, j-1, k))*ym;
        kp.values[kp.idx+3] = (j+1 == t.ny || not isSolid(mesh, i, j+1, k))*yp;
        kp.values[kp.idx+4] = (k-1 == -1   || not isSolid(mesh, i, j, k-1))*zm;
        kp.values[kp.idx+5] = (k+1 == t.nz || not isSolid(mesh, i, j, k+1))*zp;
        kp.values[kp.idx+6] = -(xp+xm+yp+ym+zp+zm);
    }

    return 0;
}

TF_Func uint32_t k_Lap_NC_FirstNb_w_solids(const tfa::KernelParams &kp)
{
    NUM_COMPONENTS(kp, 1);
    const auto &mesh = tfa::getSMesh(kp.sim);
    const auto &t = mesh.topo;
    const auto &ijk = tfa::getCellIJKFromId(mesh, kp.row);
    auto i = ijk.i;
    auto j = ijk.j;
    auto k = ijk.k;
    const double iVol = 1.0/tfa::calcCellVolume(mesh, ijk);
    const tfa::V3 areas = tfa::calcCellFacesArea(mesh, ijk);
    const double xcoef = areas[0]*iVol;
    const double ycoef = areas[1]*iVol;
    const double zcoef = areas[2]*iVol;

    if (isSolid(mesh, i, j, k))
    {
        kp.values[kp.idx+0] = 0.0;
        kp.values[kp.idx+1] = 0.0;
        kp.values[kp.idx+2] = 0.0;
        kp.values[kp.idx+3] = 0.0;
        kp.values[kp.idx+4] = 0.0;
        kp.values[kp.idx+5] = 0.0;
        kp.values[kp.idx+6] = 0.0;
    }
    else
    {
        double xm = not isSolid(mesh, i-1, j, k)? xcoef/tfa::calcNodeDX(mesh, i, i-1) : xcoef*xcoef*2.0;
        double xp = not isSolid(mesh, i+1, j, k)? xcoef/tfa::calcNodeDX(mesh, i+1, i) : xcoef*xcoef*2.0;
        double ym = not isSolid(mesh, i, j-1, k)? ycoef/tfa::calcNodeDY(mesh, j, j-1) : ycoef*ycoef*2.0;
        double yp = not isSolid(mesh, i, j+1, k)? ycoef/tfa::calcNodeDY(mesh, j+1, j) : ycoef*ycoef*2.0;
        double zm = not isSolid(mesh, i, j, k-1)? zcoef/tfa::calcNodeDZ(mesh, k, k-1) : zcoef*zcoef*2.0;
        double zp = not isSolid(mesh, i, j, k+1)? zcoef/tfa::calcNodeDZ(mesh, k+1, k) : zcoef*zcoef*2.0;
        kp.values[kp.idx+0] = (not isSolid(mesh, i-1, j, k) || i-1 == -1  )*xm;
        kp.values[kp.idx+1] = (not isSolid(mesh, i+1, j, k) || i+1 == t.nx)*xp;
        kp.values[kp.idx+2] = (not isSolid(mesh, i, j-1, k) || j-1 == -1  )*ym;
        kp.values[kp.idx+3] = (not isSolid(mesh, i, j+1, k) || j+1 == t.ny)*yp;
        kp.values[kp.idx+4] = (not isSolid(mesh, i, j, k-1) || k-1 == -1  )*zm;
        kp.values[kp.idx+5] = (not isSolid(mesh, i, j, k+1) || k+1 == t.nz)*zp;
        kp.values[kp.idx+6] = -(xp+xm+yp+ym+zp+zm);
    }

    return 0;
}

TF_Func uint32_t k_Lap_NN_FirstNb_w_solids(const tfa::KernelParams &kp)
{
    NUM_COMPONENTS(kp, 1);

    if (tfa::getNumEntries(*kp.pat, kp.row) == 2)
    {
        kp.values[kp.idx+0] = 0.0;
        kp.values[kp.idx+1] = 0.0;
        return 0;
    }

    const auto &mesh = tfa::getSMesh(kp.sim);
    const auto &t = mesh.topo;
    const auto &N = tfa::getNodeIJKFromId(mesh, kp.row);
    auto i = N.i;
    auto j = N.j;
    auto k = N.k;
    const double iVol = 1.0/tfa::calcCellVolume(mesh, N);
    const tfa::V3 areas = tfa::calcCellFacesArea(mesh, N);
    const double xcoef = areas[0]*iVol;
    const double ycoef = areas[1]*iVol;
    const double zcoef = areas[2]*iVol;

    if (isSolid(mesh, i, j, k))
    {
        kp.values[kp.idx+0] = 0.0;
        kp.values[kp.idx+1] = 0.0;
        kp.values[kp.idx+2] = 0.0;
        kp.values[kp.idx+3] = 0.0;
        kp.values[kp.idx+4] = 0.0;
        kp.values[kp.idx+5] = 0.0;
        kp.values[kp.idx+6] = 0.0;
    }
    else
    {
        double xm = (i-1 == -1   || not isSolid(mesh, i-1, j, k))? xcoef/tfa::calcNodeDX(mesh, i, i-1) : 1.0*xcoef*xcoef*2.0;
        double xp = (i+1 == t.nx || not isSolid(mesh, i+1, j, k))? xcoef/tfa::calcNodeDX(mesh, i+1, i) : 1.0*xcoef*xcoef*2.0;
        double ym = (j-1 == -1   || not isSolid(mesh, i, j-1, k))? ycoef/tfa::calcNodeDY(mesh, j, j-1) : 1.0*ycoef*ycoef*2.0;
        double yp = (j+1 == t.ny || not isSolid(mesh, i, j+1, k))? ycoef/tfa::calcNodeDY(mesh, j+1, j) : 1.0*ycoef*ycoef*2.0;
        double zm = (k-1 == -1   || not isSolid(mesh, i, j, k-1))? zcoef/tfa::calcNodeDZ(mesh, k, k-1) : 1.0*zcoef*zcoef*2.0;
        double zp = (k+1 == t.nz || not isSolid(mesh, i, j, k+1))? zcoef/tfa::calcNodeDZ(mesh, k+1, k) : 1.0*zcoef*zcoef*2.0;
        kp.values[kp.idx+0] = (i-1 == -1   || not isSolid(mesh, i-1, j, k))*xm;
        kp.values[kp.idx+1] = (i+1 == t.nx || not isSolid(mesh, i+1, j, k))*xp;
        kp.values[kp.idx+2] = (j-1 == -1   || not isSolid(mesh, i, j-1, k))*ym;
        kp.values[kp.idx+3] = (j+1 == t.ny || not isSolid(mesh, i, j+1, k))*yp;
        kp.values[kp.idx+4] = (k-1 == -1   || not isSolid(mesh, i, j, k-1))*zm;
        kp.values[kp.idx+5] = (k+1 == t.nz || not isSolid(mesh, i, j, k+1))*zp;
        kp.values[kp.idx+6] = -(xp+xm+yp+ym+zp+zm);
    }

    return 0;
}

TF_Func uint32_t k_Lap_NN_FirstNb_w_solids_neumann(const tfa::KernelParams &kp)
{
    NUM_COMPONENTS(kp, 1);

    if (tfa::getNumEntries(*kp.pat, kp.row) == 2)
    {
        kp.values[kp.idx+0] = 0.0;
        kp.values[kp.idx+1] = 0.0;
        return 0;
    }

    const auto &mesh = tfa::getSMesh(kp.sim);
    const auto &t = mesh.topo;
    const auto &N = tfa::getNodeIJKFromId(mesh, kp.row);
    auto i = N.i;
    auto j = N.j;
    auto k = N.k;
    const double iVol = 1.0/tfa::calcCellVolume(mesh, N);
    const tfa::V3 areas = tfa::calcCellFacesArea(mesh, N);
    const double xcoef = areas[0]*iVol;
    const double ycoef = areas[1]*iVol;
    const double zcoef = areas[2]*iVol;

    if (isSolid(mesh, i, j, k))
    {
        kp.values[kp.idx+0] = 0.0;
        kp.values[kp.idx+1] = 0.0;
        kp.values[kp.idx+2] = 0.0;
        kp.values[kp.idx+3] = 0.0;
        kp.values[kp.idx+4] = 0.0;
        kp.values[kp.idx+5] = 0.0;
        kp.values[kp.idx+6] = 0.0;
    }
    else
    {
        double xm = (i-1 == -1   || not isSolid(mesh, i-1, j, k))? xcoef/tfa::calcNodeDX(mesh, i, i-1) : 0.0*xcoef*xcoef*2.0;
        double xp = (i+1 == t.nx || not isSolid(mesh, i+1, j, k))? xcoef/tfa::calcNodeDX(mesh, i+1, i) : 0.0*xcoef*xcoef*2.0;
        double ym = (j-1 == -1   || not isSolid(mesh, i, j-1, k))? ycoef/tfa::calcNodeDY(mesh, j, j-1) : 0.0*ycoef*ycoef*2.0;
        double yp = (j+1 == t.ny || not isSolid(mesh, i, j+1, k))? ycoef/tfa::calcNodeDY(mesh, j+1, j) : 0.0*ycoef*ycoef*2.0;
        double zm = (k-1 == -1   || not isSolid(mesh, i, j, k-1))? zcoef/tfa::calcNodeDZ(mesh, k, k-1) : 0.0*zcoef*zcoef*2.0;
        double zp = (k+1 == t.nz || not isSolid(mesh, i, j, k+1))? zcoef/tfa::calcNodeDZ(mesh, k+1, k) : 0.0*zcoef*zcoef*2.0;
        kp.values[kp.idx+0] = (i-1 == -1   || not isSolid(mesh, i-1, j, k))*xm;
        kp.values[kp.idx+1] = (i+1 == t.nx || not isSolid(mesh, i+1, j, k))*xp;
        kp.values[kp.idx+2] = (j-1 == -1   || not isSolid(mesh, i, j-1, k))*ym;
        kp.values[kp.idx+3] = (j+1 == t.ny || not isSolid(mesh, i, j+1, k))*yp;
        kp.values[kp.idx+4] = (k-1 == -1   || not isSolid(mesh, i, j, k-1))*zm;
        kp.values[kp.idx+5] = (k+1 == t.nz || not isSolid(mesh, i, j, k+1))*zp;
        kp.values[kp.idx+6] = -(xp+xm+yp+ym+zp+zm);
    }

    return 0;
}

TF_Func uint32_t k_Interp_CF_VW_ZB_ZINT_FirstNb_w_solids(const tfa::KernelParams &k)
{
    NUM_COMPONENTS(k, 1);
    const auto &mesh = tfa::getSMesh(k.sim);
    const auto &pat = *(k.pat);
    auto numEnts = tfa::getNumEntries(pat, k.row);
    if (numEnts == 1)
    {
        k.values[k.idx] = 0.0;
        return 0;
    }

    const auto &ijk = tfa::getFaceIJKFromId(mesh, k.row);
    const auto &nbs  = tfa::getFaceNbCellsIJK(mesh, ijk);
    double V  = tfa::calcCellVolume(mesh, nbs.first);
    double OV = tfa::calcCellVolume(mesh, nbs.second);
    double v1 = (ijk.f%2)? OV : V;
    double v2 = (ijk.f%2)? V : OV;
    bool firstIsSolid = isSolid(mesh, nbs.first.i, nbs.first.j, nbs.first.k);
    // bool secondIsSolid = isSolid(mesh, nbs.second.i, nbs.second.j, nbs.second.k);
    bool secondIsSolid = tfa::isNone(nbs.second) ? firstIsSolid : isSolid(mesh, nbs.second.i, nbs.second.j, nbs.second.k);
    double w = (firstIsSolid == secondIsSolid)? 1.0 : 0.0;

    k.values[k.idx+0] = w*v1/(v1+v2);
    k.values[k.idx+1] = w*v2/(v1+v2);

    return 0;
}

TF_Func uint32_t k_Interp_NF_VW_ZB_ZINT_FirstNb_w_solids(const tfa::KernelParams &k)
{
    NUM_COMPONENTS(k, 1);

    const auto &m = tfa::getSMesh(k.sim);
    const auto &F = tfa::getFaceIJKFromId(m, k.row);

    if (tfa::isBoundaryFace(m, F))
    {
        k.values[k.idx+0] = 0.0;
        k.values[k.idx+1] = 0.0;
    }
    else
    {
        const auto &nbs = tfa::getFaceNbCellsIJK(m, F);
        const double V = tfa::calcCellVolume(m, nbs.first);
        const double OV = tfa::calcCellVolume(m, nbs.second);
        const double v1 = (F.f%2)? OV : V;
        const double v2 = (F.f%2)? V : OV;
        bool firstIsSolid = isSolid(m, nbs.first.i, nbs.first.j, nbs.first.k);
        // bool secondIsSolid = isSolid(m, nbs.second.i, nbs.second.j, nbs.second.k);
        bool secondIsSolid = tfa::isNone(nbs.second) ? firstIsSolid : isSolid(m, nbs.second.i, nbs.second.j, nbs.second.k);
        double w = (firstIsSolid == secondIsSolid)? 1.0 : 0.0;
        k.values[k.idx+0] = w*v1/(V+OV);
        k.values[k.idx+1] = w*v2/(V+OV);
    }

    return 0;
}

TF_Func uint32_t k_Interp_NF_VW_CB_ZINT_FirstNb_w_solids(const tfa::KernelParams &k)
{
    NUM_COMPONENTS(k, 1);

    const auto &m = tfa::getSMesh(k.sim);
    const auto &F = tfa::getFaceIJKFromId(m, k.row);

    if (tfa::isBoundaryFace(m, F))
    {
        k.values[k.idx+0] = 0.0;
        k.values[k.idx+1] = 1.0;
    }
    else
    {
        const auto &nbs = tfa::getFaceNbCellsIJK(m, F);
        const double V = tfa::calcCellVolume(m, nbs.first);
        const double OV = tfa::calcCellVolume(m, nbs.second);
        const double v1 = (F.f%2)? OV : V;
        const double v2 = (F.f%2)? V : OV;
        bool firstIsSolid = isSolid(m, nbs.first.i, nbs.first.j, nbs.first.k);
        // bool secondIsSolid = isSolid(m, nbs.second.i, nbs.second.j, nbs.second.k);
        bool secondIsSolid = tfa::isNone(nbs.second) ? firstIsSolid : isSolid(m, nbs.second.i, nbs.second.j, nbs.second.k);
        double w = (firstIsSolid == secondIsSolid)? 1.0 : 0.0;
        k.values[k.idx+0] = w*v1/(V+OV);
        k.values[k.idx+1] = w*v2/(V+OV);
    }

    return 0;
}

TF_Func uint32_t k_Interp_NF_SP_FirstNb_w_solids(const tfa::KernelParams &k)
{
    NUM_COMPONENTS(k, 1);
    const auto &mesh = tfa::getSMesh(k.sim);
    const auto &ijk = tfa::getFaceIJKFromId(mesh, k.row);
    const auto &nbs  = tfa::getFaceNbCellsIJK(mesh, ijk);
    bool firstIsSolid = isSolid(mesh, nbs.first.i, nbs.first.j, nbs.first.k);
    // bool secondIsSolid = isSolid(mesh, nbs.second.i, nbs.second.j, nbs.second.k);
    bool secondIsSolid = tfa::isNone(nbs.second) ? firstIsSolid : isSolid(mesh, nbs.second.i, nbs.second.j, nbs.second.k);
    double w = (firstIsSolid == secondIsSolid)? 1.0 : 0.0;
    const bool bnd = tfa::isBoundaryFace(mesh, ijk);
    k.values[k.idx+0] = (bnd)? 0.0 : w*0.5;
    k.values[k.idx+1] = (bnd)? 1.0 : w*0.5;

    return 0;
}

static uint32_t k_GradComp_CF(const tfa::KernelParams &k, int32_t comp)
{
    NUM_COMPONENTS(k, 1);
    const auto &mesh = tfa::getSMesh(k.sim);
    if (tfa::getNumEntries(*k.pat, k.row) == 1)
    {
        k.values[k.idx] = 0.0;
    }
    else
    {
        const auto &F = tfa::getFaceIJKFromId(mesh, k.row);
        const auto &nbs  = tfa::getFaceNbCellsIJK(mesh, F);
        bool firstIsSolid = isSolid(mesh, nbs.first.i, nbs.first.j, nbs.first.k);
        // bool secondIsSolid = isSolid(mesh, nbs.second.i, nbs.second.j, nbs.second.k);
        bool secondIsSolid = tfa::isNone(nbs.second) ? firstIsSolid : isSolid(mesh, nbs.second.i, nbs.second.j, nbs.second.k);
        double w = (firstIsSolid == secondIsSolid)? 1.0 : 0.0;
        double id = tfa::calcFaceArea(mesh, F)/tfa::calcFaceStaggeredVol(mesh, F);
        k.values[k.idx+0] = (F.f/2 == comp)? -id*w : 0.0;
        k.values[k.idx+1] = (F.f/2 == comp)? +id*w : 0.0;
    }
    return 0;
}

TF_Func uint32_t k_GradX_CF_FirstNb_w_solids(const tfa::KernelParams &k)
{
    return k_GradComp_CF(k, 0);
}

TF_Func uint32_t k_GradY_CF_FirstNb_w_solids(const tfa::KernelParams &k)
{
    return k_GradComp_CF(k, 1);
}

TF_Func uint32_t k_GradZ_CF_FirstNb_w_solids(const tfa::KernelParams &k)
{
    return k_GradComp_CF(k, 2);
}

static uint32_t k_GradComp_NF(const tfa::KernelParams &k, int32_t comp)
{
    NUM_COMPONENTS(k, 1);
    const auto &mesh = tfa::getSMesh(k.sim);
    const auto &F = tfa::getFaceIJKFromId(mesh, k.row);
    const auto &nbs  = tfa::getFaceNbCellsIJK(mesh, F);
    bool firstIsSolid = isSolid(mesh, nbs.first.i, nbs.first.j, nbs.first.k);
    bool secondIsSolid = tfa::isNone(nbs.second) ? firstIsSolid : isSolid(mesh, nbs.second.i, nbs.second.j, nbs.second.k);
    double w = (firstIsSolid == secondIsSolid)? 1.0 : 0.0;
    const double sign = (F.f%2==0 && tfa::isBoundaryFace(mesh, F))? -1.0: 1.0;
    double id = sign*tfa::calcFaceArea(mesh, F)/tfa::calcFaceStaggeredVol(mesh, F);
    k.values[k.idx+0] = (F.f/2 == comp)? -id*w : 0.0;
    k.values[k.idx+1] = (F.f/2 == comp)? +id*w : 0.0;
    return 0;
}

TF_Func uint32_t k_GradX_NF_FirstNb_w_solids(const tfa::KernelParams &k)
{
    return k_GradComp_NF(k, 0);
}

TF_Func uint32_t k_GradY_NF_FirstNb_w_solids(const tfa::KernelParams &k)
{
    return k_GradComp_NF(k, 1);
}

TF_Func uint32_t k_GradZ_NF_FirstNb_w_solids(const tfa::KernelParams &k)
{
    return k_GradComp_NF(k, 2);
}
