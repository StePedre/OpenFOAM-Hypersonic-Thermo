#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm> // Per std::max
#include <mutation++.h>

using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Kinetics;

// --- PARAMETRI DEL PAPER (Sezione 3.4) ---
const double T_tr_0  = 30000.0; // K
const double T_vib_0 = 1000.0;  // K
const double n_0     = 5.0e22;  // particelle/m^3 per ogni specie

// Costanti fisiche
// Nota: Usiamo Mutation::NA se disponibile o definiamo la nostra evitando conflitti
// Per sicurezza, usiamo un nome univoco per evitare 'ambiguous reference'
const double MY_NA = 6.02214076e23; 
const double MY_Ru = 8.31446;       

int main() {
    std::cout << "Simulating Chemically Reacting Air (N2-N) - Figure 7 Case" << std::endl;

    // 1. Inizializzazione Miscela
    Mutation::MixtureOptions opts("air_5");
    opts.setStateModel("ChemNonEqTTv");
    // Impostiamo il database termodinamico RRHO (Rigid Rotor Harmonic Oscillator)
    // che è standard per alte temperature e coerenza con il modello a 2 temperature.
    opts.setThermodynamicDatabase("RRHO"); 
    Mixture mix(opts);

    int ns = mix.nSpecies();
    int i_N2 = mix.speciesIndex("N2");
    int i_N  = mix.speciesIndex("N");

    // 2. Condizioni Iniziali (Paper Sezione 3.4)
    double T_tr  = T_tr_0;
    double T_vib = T_vib_0;
    
    // Densità parziali iniziali [kg/m^3]
    // n = particelle/m^3 -> rho = n * Mw / NA
    // mix.speciesMw(i) restituisce il peso molecolare in kg/mol
    std::vector<double> rho_i(ns, 0.0);
    
    // Usiamo MY_NA per evitare ambiguità con Mutation::NA
    rho_i[i_N2] = (n_0 / MY_NA) * mix.speciesMw(i_N2);
    rho_i[i_N]  = (n_0 / MY_NA) * mix.speciesMw(i_N);

    // Calcolo densità totale
    double rho = 0.0;
    for(double r : rho_i) rho += r;

    // Set stato iniziale
    // In ChemNonEqTTv, il vettore T ha 2 componenti: [T_heavy, T_internal]
    std::vector<double> T_vec = {T_tr, T_vib};
    
    // setState accetta: (densità_parziali, temperature, stato_variabili)
    // stato_variabili = 1 significa che stiamo settando (rho_i, T)
    mix.setState(rho_i.data(), T_vec.data(), 1); 

    // Calcolo Energia Totale per unità di volume (che deve conservarsi in un sistema adiabatico isocoro)
    // E = u [J/kg] * rho [kg/m^3] -> [J/m^3]
    double E_tot_conserved_density = mix.mixtureEnergyMass() * rho;

    std::cout << "Initial State:" << std::endl;
    std::cout << "  Rho_N2: " << rho_i[i_N2] << " kg/m3" << std::endl;
    std::cout << "  Rho_N : " << rho_i[i_N] << " kg/m3" << std::endl;
    std::cout << "  T_tr  : " << T_tr << " K" << std::endl;
    std::cout << "  T_vib : " << T_vib << " K" << std::endl;
    std::cout << "  E_tot : " << E_tot_conserved_density << " J/m3" << std::endl;

    // Output CSV a video per plotting
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    std::cout << "Time|T_tr|T_vib|N_N2|N_N" << std::endl;

    // Vettori di lavoro
    std::vector<double> wdot(ns);          // Ratei di produzione [kg/m^3 s]
    std::vector<double> source_en(mix.nEnergyEqns()); // Termini sorgente energia [J/m^3 s]
    
    double time = 0.0;
    // Timestep iniziale piccolo per catturare la rapida cinetica iniziale
    double dt = 1.0e-13; 
    double t_end = 1.0e-3;
    int step = 0;
    
    // Variabili temporanee per i loop di Newton
    std::vector<double> T_iter(2);
    std::vector<double> T_pert(2);

    while(time < t_end) {
        
        // --- 1. AGGIORNAMENTO STATO CORRENTE ---
        T_vec[0] = T_tr;
        T_vec[1] = T_vib;
        mix.setState(rho_i.data(), T_vec.data(), 1);
        
        // --- 2. CALCOLO RATEI ---
        // wdot = produzione netta specie i (chimica) [kg/m^3 s]
        mix.netProductionRates(wdot.data()); 
        
        // source_en = termine sorgente energia [J/m^3 s]
        // source_en[0] è la sorgente per l'equazione dell'energia vibrazionale
        // Include: Scambio V-T (Millikan-White) + Energia persa/guadagnata da reazioni (Chemistry-Vibration coupling)
        mix.energyTransferSource(source_en.data());

        // --- 3. INTEGRAZIONE DENSITA' ---
        // rho_i(t+dt) = rho_i(t) + wdot_i * dt
        for(int i=0; i<ns; ++i) {
            rho_i[i] += wdot[i] * dt;
            if(rho_i[i] < 0.0) rho_i[i] = 1.0e-20; // Protezione numeri negativi
        }
        
        // Aggiorna densità totale
        rho = 0.0;
        for(double r : rho_i) rho += r;

        // --- 4. INTEGRAZIONE TEMPERATURA VIBRAZIONALE ---
        // Calcoliamo il calore specifico vibrazionale volumetrico Cv_vib [J/m^3 K]
        // Lo facciamo perturbando T_vib e vedendo quanto cambia l'energia interna della miscela.
        // Nota: manteniamo le NUOVE densità rho_i costanti durante questa perturbazione.
        
        double deltaT = 10.0; // Perturbazione K
        
        // Stato base (T_tr, T_vib) con nuove densità
        T_iter[0] = T_tr; T_iter[1] = T_vib;
        mix.setState(rho_i.data(), T_iter.data(), 1);
        double E_base = mix.mixtureEnergyMass() * rho;
        
        // Stato perturbato (T_tr, T_vib + delta)
        T_pert[0] = T_tr; T_pert[1] = T_vib + deltaT;
        mix.setState(rho_i.data(), T_pert.data(), 1);
        double E_pert = mix.mixtureEnergyMass() * rho;
        
        double Cv_vib_vol = (E_pert - E_base) / deltaT;
        
        // Equazione energia: d(E_vib_vol)/dt = source_en[0]
        // Approssimazione: rho * Cv_vib * dT_vib/dt = source_en[0]
        // Nota: questa è un'approssimazione che trascura il termine d(rho)/dt * e_vib. 
        // Per heat bath a rho=costante è esatta, ma qui rho_i cambia.
        // Tuttavia, source_en[0] in Mutation++ di solito include già i termini di accoppiamento chimico.
        // Quindi integriamo direttamente T_vib.
        
        double dT_vib = (source_en[0] * dt) / (Cv_vib_vol + 1.0e-20);
        T_vib += dT_vib;
        
        // Clamp T_vib per stabilità
        if(T_vib < 100.0) T_vib = 100.0;
        if(T_vib > 40000.0) T_vib = 40000.0;

        // --- 5. INTEGRAZIONE TEMPERATURA TRASLAZIONALE ---
        // Usiamo la conservazione dell'energia totale per trovare T_tr.
        // E_tot_initial = E_internal(T_tr_new, T_vib_new, rho_i_new)
        // Risolviamo per T_tr_new usando Newton-Raphson.
        
        for(int k=0; k<10; ++k) {
            T_iter[0] = T_tr; T_iter[1] = T_vib;
            mix.setState(rho_i.data(), T_iter.data(), 1);
            
            double E_curr = mix.mixtureEnergyMass() * rho;
            double diff = E_curr - E_tot_conserved_density;
            
            // Convergenza?
            if(std::abs(diff) < 1.0e-3 * E_tot_conserved_density) break;
            
            // Calcolo derivata dE/dT_tr (Cv_tr volumetrico)
            T_pert[0] = T_tr + deltaT; T_pert[1] = T_vib;
            mix.setState(rho_i.data(), T_pert.data(), 1);
            double E_p = mix.mixtureEnergyMass() * rho;
            double Cv_tr_vol = (E_p - E_curr) / deltaT;
            
            // Newton step: T_new = T_old - f(T)/f'(T)
            T_tr -= diff / (Cv_tr_vol + 1.0e-20);
            
            if(T_tr < 100.0) T_tr = 100.0; // Clamp
        }

        // --- 6. ADATTAMENTO TIMESTEP ---
        // Aumentiamo il timestep man mano che il sistema rallenta
        if (time > 1e-10 && dt < 1e-11) dt = 1.0e-11;
        if (time > 1e-9 && dt < 1.0e-10) dt = 1.0e-10;
        if (time > 1e-8 && dt < 1.0e-9)  dt = 1.0e-9;
        if (time > 1e-7 && dt < 5.0e-9)  dt = 5.0e-9;
        if (time > 1e-6 && dt < 2.0e-8)  dt = 2.0e-8;
dt=1.0e-12;
        time += dt;
        step++;

        // Output ogni tot passi
        if(step % 500 == 0) {
            // Calcolo densità numeriche normalizzate
            double n_N2 = (rho_i[i_N2] / mix.speciesMw(i_N2)) * MY_NA;
            double n_N  = (rho_i[i_N]  / mix.speciesMw(i_N)) * MY_NA;
            
            std::cout.precision(6);
            std::cout << std::scientific << time << " | " 
                      << std::fixed << T_tr << " | " 
                      << T_vib << " | " 
                      << n_N2 / n_0 << " | " 
                      << n_N / n_0 
                      << std::endl;
        }
    }

    std::cout << "Simulazione completata." << std::endl;
    return 0;
}
