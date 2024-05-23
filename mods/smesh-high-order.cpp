#include "tfa/Simulation.h"

using SM = tfa::SMesh;

static void periodicWrapAround(const SM &m, int32_t &i, int32_t &j, int32_t &k)
{
    const auto &t = m.topo;
    if (isPeriodicX(m)) i = (i%t.nx)+t.nx*(i<0);
    if (isPeriodicY(m)) j = (j%t.ny)+t.ny*(j<0);
    if (isPeriodicZ(m)) k = (k%t.nz)+t.nz*(k<0);
}

static int64_t my_getCellIndex(const SM &m, int32_t i, int32_t j, int32_t k)
{
    periodicWrapAround(m, i, j, k);
    const auto &t = m.topo;
    if (i < 0 || i >= t.nx || j < 0 || j >= t.ny || k < 0 || k >= t.nz) return -1;
    auto p = getCellPartition(t, i, j, k);
    const auto i0 = t.partitionBlocks[p].i0;
    const auto j0 = t.partitionBlocks[p].j0;
    const auto k0 = t.partitionBlocks[p].k0;
    const auto px = t.partitionBlocks[p].px;
    const auto py = t.partitionBlocks[p].py;
    int64_t result = (k-k0)*px*py+(j-j0)*px+(i-i0)+TO_S64(t.firstOwnedCells[p]);
    return result;
}

TF_Func tfa::Pattern p_CC_cross(const tfa::PatternParams &p)
{
    const auto &mesh = tfa::getSMesh(p.sim);
    tfa::Pattern result;
    result.srcPtr = &tfa::getDomain(p.sim, "Cells");
    result.tgtPtr = &tfa::getDomain(p.sim, "Cells");

    uint32_t rows = tfa::getNumCells(mesh);
    uint64_t cols = tfa::getNumCells(mesh);
    auto &pat = result.pat;
    pat.offset.reserve(rows+1);
    pat.offset.push_back(0);
    pat.numCols = tfa::allReduce(cols, MPI_SUM);
    for (uint32_t c = 0; c < rows; c++)
    {
        auto ijk = tfa::getCellIJKFromId(mesh, c);
        auto i = ijk.i;
        auto j = ijk.j;
        auto k = ijk.k;

        if (my_getCellIndex(mesh, i-2, j, k) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i-2, j, k));
        if (my_getCellIndex(mesh, i-1, j, k) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i-1, j, k));
        if (my_getCellIndex(mesh, i+1, j, k) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i+1, j, k));
        if (my_getCellIndex(mesh, i+2, j, k) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i+2, j, k));

        if (my_getCellIndex(mesh, i, j-2, k) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i, j-2, k));
        if (my_getCellIndex(mesh, i, j-1, k) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i, j-1, k));
        if (my_getCellIndex(mesh, i, j+1, k) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i, j+1, k));
        if (my_getCellIndex(mesh, i, j+2, k) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i, j+2, k));

        if (my_getCellIndex(mesh, i, j, k-2) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i, j, k-2));
        if (my_getCellIndex(mesh, i, j, k-1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i, j, k-1));
        if (my_getCellIndex(mesh, i, j, k+1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i, j, k+1));
        if (my_getCellIndex(mesh, i, j, k+2) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i, j, k+2));

        pat.data.push_back(tfa::getCellIndex(mesh, i, j, k));
        pat.offset.push_back(TO_S32(pat.data.size()));
    }

    tfa::updateRequiredCols(pat, tfa::PDist(TO_LOC(cols)));
    return result;
}

TF_Func tfa::Pattern p_CC_cube(const tfa::PatternParams &p)
{
    const auto &mesh = tfa::getSMesh(p.sim);
    tfa::Pattern result;
    result.srcPtr = &tfa::getDomain(p.sim, "Cells");
    result.tgtPtr = &tfa::getDomain(p.sim, "Cells");

    uint32_t rows = tfa::getNumCells(mesh);
    uint64_t cols = tfa::getNumCells(mesh);
    auto &pat = result.pat;
    pat.offset.reserve(rows+1);
    pat.offset.push_back(0);
    pat.numCols = tfa::allReduce(cols, MPI_SUM);
    for (uint32_t c = 0; c < rows; c++)
    {
        auto ijk = tfa::getCellIJKFromId(mesh, c);
        auto i = ijk.i;
        auto j = ijk.j;
        auto k = ijk.k;

        if (my_getCellIndex(mesh, i-1, j-1, k-1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i-1, j-1, k-1));
        if (my_getCellIndex(mesh, i  , j-1, k-1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i  , j-1, k-1));
        if (my_getCellIndex(mesh, i+1, j-1, k-1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i+1, j-1, k-1));
        if (my_getCellIndex(mesh, i-1, j  , k-1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i-1, j  , k-1));
        if (my_getCellIndex(mesh, i  , j  , k-1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i  , j  , k-1));
        if (my_getCellIndex(mesh, i+1, j  , k-1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i+1, j  , k-1));
        if (my_getCellIndex(mesh, i-1, j+1, k-1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i-1, j+1, k-1));
        if (my_getCellIndex(mesh, i  , j+1, k-1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i  , j+1, k-1));
        if (my_getCellIndex(mesh, i+1, j+1, k-1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i+1, j+1, k-1));

        if (my_getCellIndex(mesh, i-1, j-1, k)   != -1) pat.data.push_back(tfa::getCellIndex(mesh, i-1, j-1, k));
        if (my_getCellIndex(mesh, i  , j-1, k)   != -1) pat.data.push_back(tfa::getCellIndex(mesh, i  , j-1, k));
        if (my_getCellIndex(mesh, i+1, j-1, k)   != -1) pat.data.push_back(tfa::getCellIndex(mesh, i+1, j-1, k));
        if (my_getCellIndex(mesh, i-1, j  , k)   != -1) pat.data.push_back(tfa::getCellIndex(mesh, i-1, j  , k));
        // NOTE: The entry for i, j, k is the last one.
        if (my_getCellIndex(mesh, i+1, j  , k)   != -1) pat.data.push_back(tfa::getCellIndex(mesh, i+1, j  , k));
        if (my_getCellIndex(mesh, i-1, j+1, k)   != -1) pat.data.push_back(tfa::getCellIndex(mesh, i-1, j+1, k));
        if (my_getCellIndex(mesh, i  , j+1, k)   != -1) pat.data.push_back(tfa::getCellIndex(mesh, i  , j+1, k));
        if (my_getCellIndex(mesh, i+1, j+1, k)   != -1) pat.data.push_back(tfa::getCellIndex(mesh, i+1, j+1, k));

        if (my_getCellIndex(mesh, i-1, j-1, k+1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i-1, j-1, k+1));
        if (my_getCellIndex(mesh, i  , j-1, k+1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i  , j-1, k+1));
        if (my_getCellIndex(mesh, i+1, j-1, k+1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i+1, j-1, k+1));
        if (my_getCellIndex(mesh, i-1, j  , k+1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i-1, j  , k+1));
        if (my_getCellIndex(mesh, i  , j  , k+1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i  , j  , k+1));
        if (my_getCellIndex(mesh, i+1, j  , k+1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i+1, j  , k+1));
        if (my_getCellIndex(mesh, i-1, j+1, k+1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i-1, j+1, k+1));
        if (my_getCellIndex(mesh, i  , j+1, k+1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i  , j+1, k+1));
        if (my_getCellIndex(mesh, i+1, j+1, k+1) != -1) pat.data.push_back(tfa::getCellIndex(mesh, i+1, j+1, k+1));

        pat.data.push_back(tfa::getCellIndex(mesh, i, j, k));
        pat.offset.push_back(TO_S32(pat.data.size()));
    }

    tfa::updateRequiredCols(pat, tfa::PDist(TO_LOC(cols)));
    return result;
}

TF_Func uint32_t k_Vol_Lap_CC_cross(const tfa::KernelParams &kp)
{
    NUM_COMPONENTS(kp, 1);
    const auto &mesh = tfa::getSMesh(kp.sim);
    const auto &ijk = tfa::getCellIJKFromId(mesh, kp.row);
    const double iVol = 1.0;
    const tfa::V3 areas = tfa::calcCellFacesArea(mesh, ijk);
    const double xcoef = areas[0]*iVol;
    const double ycoef = areas[1]*iVol;
    const double zcoef = areas[2]*iVol;
    uint32_t idx = kp.idx;
    double coef = 0.0;
    double fakeCoef = 0.0;

    auto i = ijk.i;
    auto j = ijk.j;
    auto k = ijk.k;

    if (my_getCellIndex(mesh, i-2, j, k) != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i-1, j, k) != -1) coef += (kp.values[idx++] = xcoef/tfa::calcCellDX(mesh, i,   i-1));
    if (my_getCellIndex(mesh, i+1, j, k) != -1) coef += (kp.values[idx++] = xcoef/tfa::calcCellDX(mesh, i+1, i));
    if (my_getCellIndex(mesh, i+2, j, k) != -1) coef += (kp.values[idx++] = fakeCoef);

    if (my_getCellIndex(mesh, i, j-2, k) != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i, j-1, k) != -1) coef += (kp.values[idx++] = ycoef/tfa::calcCellDY(mesh, j,   j-1));
    if (my_getCellIndex(mesh, i, j+1, k) != -1) coef += (kp.values[idx++] = ycoef/tfa::calcCellDY(mesh, j+1, j));
    if (my_getCellIndex(mesh, i, j+2, k) != -1) coef += (kp.values[idx++] = fakeCoef);

    if (my_getCellIndex(mesh, i, j, k-2) != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i, j, k-1) != -1) coef += (kp.values[idx++] = zcoef/tfa::calcCellDZ(mesh, k,   k-1));
    if (my_getCellIndex(mesh, i, j, k+1) != -1) coef += (kp.values[idx++] = zcoef/tfa::calcCellDZ(mesh, k+1, k));
    if (my_getCellIndex(mesh, i, j, k+2) != -1) coef += (kp.values[idx++] = fakeCoef);

    kp.values[idx] = -coef;

    return 0;
}

TF_Func uint32_t k_Vol_Lap_CC_cube(const tfa::KernelParams &kp)
{
    NUM_COMPONENTS(kp, 1);
    const auto &mesh = tfa::getSMesh(kp.sim);
    const auto &ijk = tfa::getCellIJKFromId(mesh, kp.row);
    const double iVol = 1.0;
    const tfa::V3 areas = tfa::calcCellFacesArea(mesh, ijk);
    const double xcoef = areas[0]*iVol;
    const double ycoef = areas[1]*iVol;
    const double zcoef = areas[2]*iVol;
    uint32_t idx = kp.idx;
    double coef = 0.0;
    double fakeCoef = 0.0;

    auto i = ijk.i;
    auto j = ijk.j;
    auto k = ijk.k;

    if (my_getCellIndex(mesh, i-1, j-1, k-1) != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i  , j-1, k-1) != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i+1, j-1, k-1) != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i-1, j  , k-1) != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i  , j  , k-1) != -1) coef += (kp.values[idx++] = zcoef/tfa::calcCellDZ(mesh, k,   k-1));
    if (my_getCellIndex(mesh, i+1, j  , k-1) != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i-1, j+1, k-1) != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i  , j+1, k-1) != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i+1, j+1, k-1) != -1) coef += (kp.values[idx++] = fakeCoef);

    if (my_getCellIndex(mesh, i-1, j-1, k)   != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i  , j-1, k)   != -1) coef += (kp.values[idx++] = ycoef/tfa::calcCellDY(mesh, j,   j-1));
    if (my_getCellIndex(mesh, i+1, j-1, k)   != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i-1, j  , k)   != -1) coef += (kp.values[idx++] = xcoef/tfa::calcCellDX(mesh, i,   i-1));
    // NOTE: The entry for i, j, k is inserted last so we have the sum of all
    // other entries.
    if (my_getCellIndex(mesh, i+1, j  , k)   != -1) coef += (kp.values[idx++] = xcoef/tfa::calcCellDX(mesh, i+1, i));
    if (my_getCellIndex(mesh, i-1, j+1, k)   != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i  , j+1, k)   != -1) coef += (kp.values[idx++] = ycoef/tfa::calcCellDY(mesh, j+1, j));
    if (my_getCellIndex(mesh, i+1, j+1, k)   != -1) coef += (kp.values[idx++] = fakeCoef);

    if (my_getCellIndex(mesh, i-1, j-1, k+1) != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i  , j-1, k+1) != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i+1, j-1, k+1) != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i-1, j  , k+1) != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i  , j  , k+1) != -1) coef += (kp.values[idx++] = zcoef/tfa::calcCellDZ(mesh, k+1, k));
    if (my_getCellIndex(mesh, i+1, j  , k+1) != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i-1, j+1, k+1) != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i  , j+1, k+1) != -1) coef += (kp.values[idx++] = fakeCoef);
    if (my_getCellIndex(mesh, i+1, j+1, k+1) != -1) coef += (kp.values[idx++] = fakeCoef);

    kp.values[idx] = -coef;

    return 0;
}
