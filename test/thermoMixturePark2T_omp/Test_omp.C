/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Example of coupling between OpenFOAM and Mutation++
    Case: Non-Reacting N2 Relaxation (Zero-Dimensional Heat Bath)
    OpenMP: run many independent cells in parallel
    Reference: aerospace-03-00034, Figure 3a
\*---------------------------------------------------------------------------*/

#include "mutation++.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <omp.h>

#include "IFstream.H"
#include "ODESolver.H"
#include "ODESystem.H"
#include "dictionary.H"

static inline std::string makeFileName(const int cellId)
{
    std::ostringstream oss;
    oss << "results_heatbath_N2_cell" << std::setw(4) << std::setfill('0') << cellId << ".dat";
    return oss.str();
}

struct MixConst
{
    int i_N2;
    int i_N;
    double m_N2;
    double inv_m_N2;
    double inv_m_N2_half;
};

static inline MixConst buildMixConst(Mutation::Mixture &mix)
{
    MixConst mc{};
    mc.i_N2 = mix.speciesIndex("N2");
    mc.i_N = mix.speciesIndex("N");
    mc.m_N2 = mix.speciesMw(mc.i_N2) / Mutation::NA;
    mc.inv_m_N2 = 1.0 / mc.m_N2;
    mc.inv_m_N2_half = 1.0 / (mc.m_N2 * 0.5);
    return mc;
}

void test2Cell(const int cellId,
               const bool writeFile,
               Mutation::Mixture &mix,
               const MixConst &mc)
{
    // Initial conditions
    double T_tr = 30000.0;
    double T_v = 1000.0;
    const double P_init = 1.0 * Mutation::ONEATM;

    // Constants
    const double KB = Mutation::KB;
    const double RU = 8.314462618;
    const double theta_v = 3371.0;
    const double D_N2_J = 9.7593 * 1.60218e-19;
    const double inv_ONEATM = 1.0 / Mutation::ONEATM;

    // Densities
    const double n_total_init = P_init / (KB * T_tr);
    double rho_N2 = n_total_init * mc.m_N2;
    double rho_N = 0.0;

    double dt = 1.0e-13;
    double t_curr = 0.0;

    const double print_interval = 1.0e-9;
    double next_print_time = 0.0;

    std::ofstream file;
    if (writeFile)
    {
        file.open(makeFileName(cellId));
        file << "Time(s) T_tr(K) T_ve(K) rho_N2 rho_N\n";
        file << 0.0 << " " << T_tr << " " << T_v << " " << rho_N2 << " " << rho_N << "\n";
    }

    // Work arrays (thread-local)
    double wdot[5] = {0.0};
    double rho_i[5] = {0.0};
    double T_vec[2] = {0.0, 0.0};
    double ev_N2 = 0.0;

    // Precomputations
    const double A_mw = 1.16e-3 * std::sqrt(mc.m_N2 * Mutation::NA) * std::cbrt(theta_v) * theta_v;
    const double B_mw = 0.015 * std::pow(mc.m_N2 * Mutation::NA, 0.25);

    const double R_N2 = RU / (mc.m_N2 * Mutation::NA);
    const double R_N = RU / (mc.m_N2 * 0.5 * Mutation::NA);

    const double v_bar_coeff = std::sqrt(8.0 * KB / (3.14159 * mc.m_N2));

    while (t_curr < 1.0e-3)
    {
        // Mutation++ state + chemistry
        rho_i[mc.i_N2] = rho_N2;
        rho_i[mc.i_N] = rho_N;
        T_vec[0] = T_tr;
        T_vec[1] = T_v;

        mix.setState(rho_i, T_vec, 1);
        mix.netProductionRates(wdot);

        // Properties
        const double n_N2 = rho_N2 * mc.inv_m_N2;
        const double n_N = rho_N * mc.inv_m_N2_half;
        const double n_tot = n_N2 + n_N;
        const double P_curr = n_tot * KB * T_tr;

        // Vibrational energy
        const double denom = std::exp(theta_v / T_v) - 1.0;
        ev_N2 = KB * theta_v / denom;
        const double ev_eq = KB * theta_v / (std::exp(theta_v / T_tr) - 1.0);

        // Chemistry rates
        const double rate_N2_kg = wdot[mc.i_N2];
        const double rate_N2_part = rate_N2_kg * mc.inv_m_N2;

        // V-T relaxation
        const double P_atm = P_curr * inv_ONEATM;
        const double arg = A_mw * (1.0 / std::cbrt(T_tr) - B_mw) - 18.42;
        const double tau_MW = (1.0 / P_atm) * std::exp(arg);

        const double temp_ratio = 50000.0 / T_tr;
        const double sigma = 1e-21 * (temp_ratio * temp_ratio);

        const double v_bar = v_bar_coeff * std::sqrt(T_tr);
        const double tau_P = 1.0 / (sigma * v_bar * n_tot);
        const double tau = tau_MW + tau_P;
        const double inv_tau = 1.0 / tau;

        // Semi-implicit integration
        const double max_change = 0.1 * rho_N2;
        const double abs_rate_N2 = std::abs(rate_N2_kg);
        if (abs_rate_N2 * dt > max_change)
            dt = max_change / abs_rate_N2;

        rho_N2 += rate_N2_kg * dt;
        rho_N += wdot[mc.i_N] * dt;

        if (rho_N2 <= 1.0e-30)
            rho_N2 = 1.0e-30;

        // Vibration source terms
        const double E_rem = 0.3 * D_N2_J;
        const double Q_VT = n_N2 * (ev_eq - ev_N2) * inv_tau;
        const double Q_CV = rate_N2_part * (E_rem - ev_N2);

        double d_ev = (Q_VT + Q_CV) / (n_N2 + 1.0e-30) * dt;

        const double limit_ev = 0.2 * ev_N2;
        if (std::abs(d_ev) > limit_ev)
            d_ev = (d_ev > 0 ? limit_ev : -limit_ev);

        ev_N2 += d_ev;
        if (ev_N2 < 1e-25)
            ev_N2 = 1e-25;

        T_v = theta_v / std::log((KB * theta_v) / ev_N2 + 1.0);

        // Translational update
        const double Cv_vol_mix = rho_N2 * 2.5 * R_N2 + rho_N * 1.5 * R_N;
        const double E_rem_tr = D_N2_J - E_rem;
        const double Chem_Power_Tr = rate_N2_part * E_rem_tr;

        const double d_Ttr_dt = (-Q_VT + Chem_Power_Tr) / Cv_vol_mix;
        T_tr += d_Ttr_dt * dt;

        // Log
        if (writeFile && t_curr >= next_print_time)
        {
            file << t_curr << " " << T_tr << " " << T_v << " " << rho_N2 << " " << rho_N << "\n";
            next_print_time += print_interval;
        }

        // timestep adaptation
        const double rate_T = std::abs(d_Ttr_dt * dt) / T_tr;
        if (rate_T < 0.001 && dt < 1e-8)
            dt *= 1.1;
        if (dt < 1e-14)
            dt = 1e-14;

        t_curr += dt;
    }
}

int main(int argc, char *argv[])
{
    const int nCells = (argc > 1) ? std::atoi(argv[1]) : 64;
    const bool writeFiles = (argc > 2) ? (std::atoi(argv[2]) != 0) : false;

    // Decide how many threads will be used (from OMP_NUM_THREADS or runtime default)
    const int nThreads = omp_get_max_threads();

    // Pre-create one Mixture per thread IN SERIAL (avoid ctor/dtor races)
    std::vector<Mutation::Mixture *> mixes(nThreads, nullptr);
    std::vector<MixConst> mconst(nThreads);

    for (int t = 0; t < nThreads; ++t)
    {
        Mutation::MixtureOptions opts("air_5");
        opts.setStateModel("ChemNonEqTTv");
        mixes[t] = new Mutation::Mixture(opts); // intentionally leaked (no delete)
        mconst[t] = buildMixConst(*mixes[t]);
    }

    const double t0 = omp_get_wtime();

#pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        Mutation::Mixture &mix = *mixes[tid];
        const MixConst &mc = mconst[tid];

#pragma omp for schedule(static)
        for (int c = 0; c < nCells; ++c)
        {
            test2Cell(c, writeFiles, mix, mc);
        }
    }

    const double t1 = omp_get_wtime();

    Foam::Info << "Done. nCells=" << nCells
               << " nThreads=" << nThreads
               << " wallTime=" << (t1 - t0) << " s" << Foam::endl;

    // DO NOT delete mixes[t];  (workaround for double-free in Mutation++ teardown)
    return 0;
}
