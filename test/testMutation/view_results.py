import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import time

# --- CONFIGURAZIONE ---
LIVE_MODE = False       # <--- Imposta a True per il live plotting, False per plot statico
REFRESH_RATE = 1.0     # Secondi di attesa tra un aggiornamento e l'altro (solo se live)
FILENAME_SIM = 'results_heatbath_N2.dat'

# --- 1. CONFIGURAZIONE STILE (Academic Look) ---
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset'] = 'cm' 
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 1.0     
plt.rcParams['xtick.direction'] = 'in'   
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top'] = True         
plt.rcParams['ytick.right'] = True       

def update_plot(ax):
    """
    Funzione che legge i dati, pulisce l'asse e ridisegna il grafico.
    Restituisce True se il plot è avvenuto con successo, False se ci sono problemi col file.
    """
    # --- 2. CARICAMENTO DATI (Con gestione errori per lettura live) ---
    try:
        # Usiamo engine='python' o gestiamo eccezioni nel caso il file sia in scrittura
        df_sim = pd.read_csv(FILENAME_SIM, sep=r'\s+')
        
        # Controllo se il dataframe è vuoto
        if df_sim.empty:
            return False
            
    except (FileNotFoundError, pd.errors.EmptyDataError):
        # Se il file non esiste ancora o è vuoto (inizio simulazione), non crashare
        print(f"In attesa del file {FILENAME_SIM}...")
        return False
    except Exception as e:
        print(f"Errore lettura file: {e}")
        return False

    # Preparazione dati Simulazione
    time_sim = df_sim['Time(s)']
    T_tr_sim = df_sim['T_tr(K)'] / 1000.0 
    T_ve_sim = df_sim['T_ve(K)'] / 1000.0

    # --- 3. PULIZIA E DISEGNO ---
    ax.clear() # Pulisce il grafico precedente

    # --- A. Our solver DATA ---
    ax.plot(time_sim.values, T_tr_sim.values, color='black', linestyle='-', linewidth=1.5, label=r'Our solver')
    ax.plot(time_sim.values, T_ve_sim.values, color='red', linestyle='-', linewidth=1.5)

    # --- 4. FORMATTAZIONE ASSI (Va riapplicata dopo ax.clear) ---
    ax.set_xscale('log')
    ax.set_xlim(1e-9, 1e-3)
    ax.set_ylim(0, 31)

    ax.set_xlabel(r'Time, $t$ [s]', fontsize=14)
    ax.set_ylabel(r'Temperature [K $\times 10^3$]', fontsize=14)

    ax.minorticks_on()
    ax.tick_params(axis='both', which='major', length=6, width=1)
    ax.tick_params(axis='both', which='minor', length=3, width=0.5)

    # --- 5. LEGENDA ---
    legend = ax.legend(loc='center left', bbox_to_anchor=(0.05, 0.45), 
                       frameon=True, edgecolor='black', fancybox=False, fontsize=11)
    legend.get_frame().set_linewidth(1.0)
    legend.set_title(r"Reacting N$_2 + N$")
    plt.setp(legend.get_title(), fontsize=12)
    
    return True

def main():
    fig, ax = plt.subplots(figsize=(8, 6))

    if LIVE_MODE:
        print("--- Modalità LIVE PLOTTING Attiva ---")
        print("Premi Ctrl+C nel terminale per interrompere.")
        plt.ion()  # Attiva modalità interattiva
        
        while True:
            try:
                success = update_plot(ax)
                if success:
                    plt.draw()
                    plt.pause(REFRESH_RATE) # Pausa per permettere l'aggiornamento grafico
                else:
                    # Se fallisce la lettura (es. file non ancora creato), aspetta comunque
                    time.sleep(REFRESH_RATE)
            except KeyboardInterrupt:
                print("\nInterruzione manuale.")
                break
        
        plt.ioff() # Disattiva modalità interattiva alla fine
        print("Salvataggio grafico finale...")
        plt.tight_layout()
        plt.savefig('figure_3a.png', dpi=300)
        plt.show() # Mostra l'ultimo stato fermo
        
    else:
        # --- MODALITÀ STATICA (Classica) ---
        print("--- Modalità STATICA ---")
        update_plot(ax)
        plt.tight_layout()
        plt.savefig('figure_3a.png', dpi=300)
        plt.show()

if __name__ == "__main__":
    main()