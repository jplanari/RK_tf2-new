void printMaxVals(std::string point, tfa::Field &ux, tfa::Field &uy, tfa::Field &uz)
{
    tfa::info("max values at this point: %s\n",point.c_str());
    tfa::info("\tux=%e\n",tfa::max(tfa::oper_max(ux)));
    tfa::info("\tuy=%e\n",tfa::max(tfa::oper_max(uy)));
    tfa::info("\tuz=%e\n",tfa::max(tfa::oper_max(uz)));
}
/*
static void calcStaggeredVolumes(const tfa::UMesh &mesh, double *buffer)
{
    for (auto it : mesh.getFacesRange())
    {
        const auto &nbs = mesh.getFaceNbCells(it);
        const auto &fnv = mesh.getFaceNormal(it);
        tfa::V3 d;
        if (nbs.second != mt::CellIdT::None)
        {
            d = mesh.getCellCentroid(nbs.second)-mesh.getCellCentroid(nbs.first);
        }
        else
        {
            d = mesh.getFaceCentroid(it)-mesh.getCellCentroid(nbs.first);
        }
        *buffer++ = (mesh.getFaceArea(it)*fabs(fnv*d));
    }
}
*/
static void calcVolumes(const tfa::SMesh &mesh, double *buffer)
{
    const auto nc = tfa::getNumCells(mesh);
    for (uint32_t it = 0; it < nc; it++)
    {
        const auto C = tfa::getCellIJKFromId(mesh, it);
        *buffer++ = tfa::calcCellVolume(mesh, C);
    }
}

static void calcVolumes(const tfa::UMesh &mesh, double *buffer)
{
    for (auto it : mesh.getCellsRange())
    {
        *buffer++ = mesh.getCellVolume(it);
    }
}
/*
static void calcInvVolumes(const tfa::UMesh &mesh, double *buffer)
{
    for (auto it : mesh.getCellsRange())
    {
        *buffer++ = 1.0/mesh.getCellVolume(it);
    }
}
*/
static void genVolumesField(tfa::Simulation &sim)
{
    if (not tfa::hasField(sim, "Omega_C"))
    {
        std::vector<double> buffer(tfa::getDomainSize(sim, "Cells"));
        if (hasSMesh(sim))
        {
            const auto &m = tfa::getSMesh(sim);
            calcVolumes(m, buffer.data());
        }
        else
        {
            const auto &m = tfa::getUMesh(sim);
            calcVolumes(m, buffer.data());
        }
        auto &f = tfa::newField(sim, 1, "Omega_C", "Cells");
        tfa::oper_setData(f, buffer.data());
    }
}

TF_Func void SetUp_Momentum_Collocated(tfa::Simulation &sim)
{
    TF_uAssert((sim.meshes.size() == 1) xor tfa::hasSMesh(sim), "not yet implemented for multiple meshes");

    const auto &meshName = tfa::getMeshName(sim);
    const auto &msuf     = tfa::meshSuffix(meshName);
    const auto &cells    = "Cells"+msuf;
    const auto &faces    = "Faces"+msuf;
    const auto &nodes    = "Nodes"+msuf;

    const auto smesh = tfa::hasSMesh(sim);
    const auto &op   = tfa::getOfficialPath();
    const auto &pMod = op+(smesh? "/smeshPatterns.so" : "/patterns.so");
    const auto &kMod = op+(smesh? "/smeshInterpolators.so" : "/interpolators.so");
    if (not tfa::hasMatrix(sim, "Id_NC"))
    {
        tfa::newMatrix(sim, "Id_NC", "pId_NC", pMod, "kId_NC", kMod, meshName);
    }
    if (not tfa::hasMatrix(sim, "Id_CN"))
    {
        tfa::newMatrix(sim, "Id_CN", "pId_CN", pMod, "kId_CN", kMod, meshName);
    }

    genVolumesField(sim);

    tfa::getOrCreateField(sim, 1, "upredx_C",  cells);
    tfa::getOrCreateField(sim, 1, "upredy_C",  cells);
    tfa::getOrCreateField(sim, 1, "upredz_C",  cells);
    tfa::getOrCreateField(sim, 1, "momSrcx_C", cells);
    tfa::getOrCreateField(sim, 1, "momSrcy_C", cells);
    tfa::getOrCreateField(sim, 1, "momSrcz_C", cells);
    tfa::getOrCreateField(sim, 1, "P_C",       cells);

    tfa::getOrCreateField(sim, 1, "ux_F", faces);
    tfa::getOrCreateField(sim, 1, "uy_F", faces);
    tfa::getOrCreateField(sim, 1, "uz_F", faces);

 }


