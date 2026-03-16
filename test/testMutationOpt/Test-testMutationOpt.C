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

void test2Opt() {
  // Setup Mutation++
  Mutation::MixtureOptions opts("air_5");
  opts.setStateModel("ChemNonEqTTv");
  Mutation::Mixture mix(opts);

  // Pre-computation of different variables
  // - Indexes
  const int i_N2 = mix.speciesIndex("N2");
  const int i_N = mix.speciesIndex("N");
  
  // - Mass and its inverses
  const double m_N2 = mix.speciesMw(i_N2) / Mutation::NA;
  const double inv_m_N2 = 1.0 / m_N2;
  const double inv_m_N2_half = 1.0 / (m_N2 * 0.5);

  // - Constants
  const double KB = Mutation::KB;
  const double RU = 8.314462618;
  const double theta_v = 3371.0;
  const double D_N2_J = 9.7593 * 1.60218e-19;
  const double inv_ONEATM = 1.0 / Mutation::ONEATM;

  // Initial conditions
  double T_tr = 30000.0;
  double T_v = 1000.0;
  double P_init = 1.0 * Mutation::ONEATM;

  // Initial densities
  double n_total_init = P_init / (KB * T_tr);
  double rho_N2 = n_total_init * m_N2;
  double rho_N = 0.0;

  // Time variables
  double dt = 1.0e-13;
  double t_curr = 0.0;
  
  // Optimization: print logic relies on double comparison instead of costly int casts
  const double print_interval = 1.0e-9;
  double next_print_time = 0.0;

  std::ofstream file("results_heatbath_N2.dat");
  // Minor optimization: use of \n instead of std::endl to avoid flush overhead
  file << "Time(s) T_tr(K) T_ve(K) rho_N2 rho_N\n";
  file << 0.0 << " " << T_tr << " " << T_v << " " << rho_N2 << " " << rho_N << "\n";

  // Arrays for Mutation++ calls
  double wdot[5];
  double rho_i[5] = {0.0}; // Optimization: direct initialization
  double T_vec[2];
  double ev_N2 = 0.0;

  // Optimization: pre-computation of constants not changing inside the loop 
  // - Millikan-white coefficients
  const double A_mw = 1.16e-3 * std::sqrt(m_N2 * Mutation::NA) * std::cbrt(theta_v) * theta_v;  
  const double B_mw = 0.015 * std::pow(m_N2 * Mutation::NA, 0.25);

  // - Gas constant per specie
  const double R_N2 = RU / (m_N2 * Mutation::NA);
  const double R_N = RU / (m_N2 * 0.5 * Mutation::NA);

  // - Thermal velocity coefficient
  const double v_bar_coeff = std::sqrt(8.0 * KB / (3.14159 * m_N2));

  // --- MAIN LOOP ---
  while (t_curr < 1.0e-3) {
    // Update state
    rho_i[i_N2] = rho_N2;
    rho_i[i_N] = rho_N;
    T_vec[0] = T_tr;
    T_vec[1] = T_v;

    // Mutation++ state update (probable bottleneck)
    mix.setState(rho_i, T_vec, 1);

    // Properties
    // Optimization: use of multiplications instead of divisions
    double n_N2 = rho_N2 * inv_m_N2;
    double n_N = rho_N * inv_m_N2_half;
    double n_tot = n_N2 + n_N;
    double P_curr = n_tot * KB * T_tr;

    // Vibrational energy
    double denom = std::exp(theta_v / T_v) - 1.0;
    ev_N2 = KB * theta_v / denom;
    double ev_eq = KB * theta_v / (std::exp(theta_v / T_tr) - 1.0);

    // Chemistry
    mix.netProductionRates(wdot); // kg/m^3/s
    double rate_N2_kg = wdot[i_N2];
    // Optimization: multiplication by inverse
    double rate_N2_part = rate_N2_kg * inv_m_N2; 

    // D-V Relaxation
    double P_atm = P_curr * inv_ONEATM;
    // Optimization: pow(T, -1/3) -> 1.0 / cbrt(T)
    double arg = A_mw * (1.0 / std::cbrt(T_tr) - B_mw) - 18.42;
    double tau_MW = (1.0 / P_atm) * std::exp(arg);

    // Optimization: pow(x, 2) -> x*x
    double temp_ratio = 50000.0 / T_tr;
    double sigma = 1e-21 * (temp_ratio * temp_ratio);

    // Optimization: thermal velocity calculation with pre-computed coefficient
    double v_bar = v_bar_coeff * std::sqrt(T_tr);
    
    double tau_P = 1.0 / (sigma * v_bar * n_tot);
    double tau = tau_MW + tau_P;
    // Pre-compute inverse of tau for later use
    double inv_tau = 1.0 / tau; 

    // Semi-implicit integration
    double max_change = 0.1 * rho_N2;
    // Optimization: compute absolute rate once
    double abs_rate_N2 = std::abs(rate_N2_kg);

    if (abs_rate_N2 * dt > max_change)
      dt = max_change / abs_rate_N2;

    rho_N2 += rate_N2_kg * dt;
    rho_N += wdot[i_N] * dt;

    // Physical protection
    if (rho_N2 <= 1.0e-30) rho_N2 = 1.0e-30;

    // Vibrational energy update
    double E_rem = 0.3 * D_N2_J;
    // Optimization: combine Q_VT and Q_CV calculations to reduce divisions
    double Q_VT = n_N2 * (ev_eq - ev_N2) * inv_tau;
    double Q_CV = rate_N2_part * (E_rem - ev_N2);

    double d_ev = (Q_VT + Q_CV) / (n_N2 + 1.0e-30) * dt;

    // Limit the change in ev_N2 to 20% per timestep for stability
    double abs_d_ev = std::abs(d_ev);
    double limit_ev = 0.2 * ev_N2;
    
    if (abs_d_ev > limit_ev) {
      d_ev = (d_ev > 0 ? limit_ev : -limit_ev);
    }

    ev_N2 += d_ev;

    if (ev_N2 < 1e-25) ev_N2 = 1e-25;

    // Update T_v from ev_N2
    T_v = theta_v / std::log((KB * theta_v) / ev_N2 + 1.0);

    // Translational temperature update
    double Cv_vol_mix = rho_N2 * 2.5 * R_N2 + rho_N * 1.5 * R_N;
    double E_rem_tr = D_N2_J - E_rem;
    double Chem_Power_Tr = rate_N2_part * E_rem_tr;

    double d_Ttr_dt = (-Q_VT + Chem_Power_Tr) / Cv_vol_mix;
    T_tr += d_Ttr_dt * dt;

    // Logging results
    // Optimization: print logic based on double comparison to avoid int casts
    if (t_curr >= next_print_time) {
      file << t_curr << " " << T_tr << " " << T_v << " " << rho_N2 << " "
           << rho_N << "\n"; // Usa \n
      next_print_time += print_interval;
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
  test2Opt();
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  Foam::Info << "Done. Time: " << elapsed.count() << " s" << Foam::endl;
  return 0;
}