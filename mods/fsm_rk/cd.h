#include <functional>
#include <map>

void diffusive(tfa::Field &u, tfa::Field &diff, tfa::Simulation &sim)
{
    auto &L = tfa::getMatrix(sim, "Lap_NN");

    TF_uAssert(sim.IOParamD.count("kinVisc") > 0, "Parameter 'kinVisc' not defined");
    double visc = sim.IOParamD["kinVisc"];

    // Diffusive_i = visc*lap(u_i)
    tfa::oper_prod(L, u, diff, visc);
}

void diffusive(tfa::Field &u, tfa::Field &diff, double visc, tfa::Simulation &sim)
{
    auto &L = tfa::getMatrix(sim, "Lap_NN");

    // Diffusive_i = visc*lap(u_i)
    tfa::oper_prod(L, u, diff, visc);
}

void diffusiveT(tfa::Field &u, tfa::Field &diff, tfa::Simulation &sim)
{
    auto &L = tfa::getMatrix(sim, "Lap_NN");

    TF_uAssert(sim.IOParamD.count("kinVisc") > 0, "Parameter 'kinVisc' not defined");
    double visc = sim.IOParamD["kinVisc"];

    // Diffusive_i = visc*lap(u_i)
    tfa::oper_prod(L, u, diff, visc);
}

void convective(tfa::Field &u, tfa::Field &conv, tfa::Simulation &sim)
{
    auto &DX  = tfa::getMatrix(sim, "DivX_FN");
    auto &DY  = tfa::getMatrix(sim, "DivY_FN");
    auto &DZ  = tfa::getMatrix(sim, "DivZ_FN");

    auto &ufx   = tfa::getField(sim, "ux_F");
    auto &ufy   = tfa::getField(sim, "uy_F");
    auto &ufz   = tfa::getField(sim, "uz_F");

    // Temporary variables.
    auto aux = tfa::getTmpFieldS(sim, ufx);
    // Convective_i = div(u*u_i)

    
    tfa::oper_axpy(ufx,*aux,1.0,0.0);
    tfa::oper_prod(u,*aux,1.0);
    tfa::oper_prod(DX, *aux, conv);
    tfa::oper_axpy(ufy,*aux,1.0,0.0);
    tfa::oper_prod(u,*aux,1.0);
    tfa::oper_prod(DY, *aux, conv, 1.0, 1.0);
    tfa::oper_axpy(ufz,*aux,1.0,0.0);
    tfa::oper_prod(u,*aux,1.0);
    tfa::oper_prod(DZ, *aux, conv, 1.0, 1.0);
}

