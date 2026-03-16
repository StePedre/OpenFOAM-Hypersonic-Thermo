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
    Reference: aerospace-03-00034, Figure 3a
\*---------------------------------------------------------------------------*/

#include "mutation++.h"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "IFstream.H"
#include "ODESolver.H"
#include "ODESystem.H"
#include "dictionary.H"

void test2() {
  // Setup Mutation++
  Mutation::MixtureOptions opts("air_5");
  opts.setStateModel("ChemNonEqTTv");
  Mutation::Mixture mix(opts);

  // Pre-computation of different variables
  // - Indexes and mass
  int i_N2 = mix.speciesIndex("N2");
  int i_N = mix.speciesIndex("N");
  double m_N2 = mix.speciesMw(i_N2) / Mutation::NA;

  // Initial conditions
  double T_tr = 30000.0;
  double T_v = 1000.0;
  double P_init = 1.0 * Mutation::ONEATM;

  // Constants
  double KB = Mutation::KB;
  double RU = 8.314462618;
  double theta_v = 3371.0;
  double D_N2_J = 9.7593 * 1.60218e-19;

  // Initial densities
  double n_total_init = P_init / (KB * T_tr);
  double rho_N2 = n_total_init * m_N2;
  double rho_N = 0.0;

  // Time variables
  double dt = 1.0e-13;
  double t_curr = 0.0;
  double print_interval = 1.0e-9;

  // Log setup
  std::ofstream file("results_heatbath_N2.dat");
  file << "Time(s) T_tr(K) T_ve(K) rho_N2 rho_N" << std::endl;
  file << 0.0 << " " << T_tr << " " << T_v << " " << rho_N2 << " " << rho_N
       << std::endl;

  // Arrays for Mutation++ calls
  double wdot[5];
  double rho_i[5];
  double T_vec[2];
  double ev_N2 = 0.0;

  // Initialize rho_i
  for (int i = 0; i < 5; ++i)
    rho_i[i] = 0.0;

  // --- MAIN LOOP ---
  while (t_curr < 1.0e-3) {
    // Update state
    rho_i[i_N2] = rho_N2;
    rho_i[i_N] = rho_N;
    T_vec[0] = T_tr;
    T_vec[1] = T_v;
    mix.setState(rho_i, T_vec, 1);

    // Properties
    double n_N2 = rho_N2 / m_N2;
    double n_N = rho_N / (m_N2 * 0.5);
    double n_tot = n_N2 + n_N;
    double P_curr = n_tot * KB * T_tr;

    // Vibrational energy
    double denom = std::exp(theta_v / T_v) - 1.0;
    ev_N2 = KB * theta_v / denom;
    double ev_eq = KB * theta_v / (std::exp(theta_v / T_tr) - 1.0);

    // Chemistry
    mix.netProductionRates(wdot); // kg/m^3/s
    double rate_N2_kg = wdot[i_N2];
    double rate_N2_part = rate_N2_kg / m_N2;

    // D-V Relaxation
    double P_atm = P_curr / Mutation::ONEATM;
    double A_mw =
        1.16e-3 * std::sqrt(m_N2 * Mutation::NA) * std::pow(theta_v, 4.0 / 3.0);
    double B_mw = 0.015 * std::pow(m_N2 * Mutation::NA, 0.25);
    double arg = A_mw * (std::pow(T_tr, -1.0 / 3.0) - B_mw) - 18.42;
    double tau_MW = (1.0 / P_atm) * std::exp(arg);
    double sigma = 1e-21 * std::pow(50000.0 / T_tr, 2.0);
    double v_bar = std::sqrt(8.0 * KB * T_tr / (3.14159 * m_N2));
    double tau_P = 1.0 / (sigma * v_bar * n_tot);
    double tau = tau_MW + tau_P;

    // Semi-implicit integration
    // Reduce dt only for this step for chemical stability
    double max_change = 0.1 * rho_N2;
    if (std::abs(rate_N2_kg * dt) > max_change) {
      dt = max_change / std::abs(rate_N2_kg);
    }

    rho_N2 += rate_N2_kg * dt;
    rho_N += wdot[i_N] * dt;

    // Physical protection
    if (rho_N2 <= 0.0)
      rho_N2 = 1.0e-30;

    // Vibrational energy update
    double E_rem = 0.3 * D_N2_J;
    double Source_VT = ev_eq / tau;
    double Decay_VT = 1.0 / tau;
    double Q_VT = n_N2 * (ev_eq - ev_N2) / tau;
    double Q_CV = rate_N2_part * (E_rem - ev_N2);

    double d_ev = (Q_VT + Q_CV) / n_N2 * dt;

    // Limit the change in ev_N2 to 20% per timestep for stability
    if (std::abs(d_ev) > 0.2 * ev_N2) {
      d_ev = (d_ev > 0 ? 0.2 : -0.2) * ev_N2;
    }

    ev_N2 += d_ev;

    // Physical protection on ev_N2
    if (ev_N2 < 1e-25)
      ev_N2 = 1e-25;

    // Update T_v from ev_N2
    T_v = theta_v / std::log(KB * theta_v / ev_N2 + 1.0);

    // Translational temperature update
    double R_N2 = RU / (m_N2 * Mutation::NA);
    double R_N = RU / (m_N2 * 0.5 * Mutation::NA);
    double Cv_vol_mix = rho_N2 * 2.5 * R_N2 + rho_N * 1.5 * R_N;

    double E_rem_tr = D_N2_J - E_rem;
    double Chem_Power_Tr = rate_N2_part * E_rem_tr;

    double d_Ttr_dt = (-Q_VT + Chem_Power_Tr) / Cv_vol_mix;
    T_tr += d_Ttr_dt * dt;

    // Logging results
    if ((int)(t_curr / print_interval) !=
        (int)((t_curr + dt) / print_interval)) {
      file << t_curr << " " << T_tr << " " << T_v << " " << rho_N2 << " "
           << rho_N << std::endl;
    }

    double rate_T = std::abs(d_Ttr_dt * dt) / T_tr;

    if (rate_T < 0.001 && dt < 1e-8)
      dt *= 1.1;
    if (dt < 1e-14)
      dt = 1e-14;

    t_curr += dt;
  }
}

int main(int argc, char *argv[]) {
  auto start = std::chrono::high_resolution_clock::now();
  test2();
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  Foam::Info << "Done. Time: " << elapsed.count() << " s" << Foam::endl;
  return 0;
}