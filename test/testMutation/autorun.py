import subprocess
import re
import sys

# --- CONFIGURAZIONE ---
N = 50  # Numero di esecuzioni
# Inserisci qui il comando esatto. 
# Se è uno script nella cartella corrente, usa "./Test-testMutation"
COMMAND = "Test-testMutation" 
# ----------------------

times = []

print(f"Avvio di {N} esecuzioni per il comando: '{COMMAND}'\n")

for i in range(1, N + 1):
    print(f"Esecuzione {i}/{N}...", end=" ", flush=True)

    try:
        # Esegue il comando su shell. 
        # capture_output=True cattura stdout e stderr.
        # text=True decodifica i byte in stringa automaticamente.
        result = subprocess.run(COMMAND, shell=True, capture_output=True, text=True)
        
        # Unisci stdout e stderr per cercare il tempo ovunque appaia
        full_output = result.stdout + result.stderr

        # Regex: Cerca "Time:" seguito da spazi e poi un numero decimale
        match = re.search(r"Time:\s*([\d\.]+)", full_output)

        if match:
            # Estrae il gruppo catturato (il numero) e converte in float
            time_val = float(match.group(1))
            times.append(time_val)
            print(f"Fatto ({time_val} s)")
        else:
            print("\n[!] Output non riconosciuto o tempo mancante.")
            # Opzionale: stampa l'output per debug
            # print(f"Output grezzo: {full_output.strip()}")

    except Exception as e:
        print(f"\n[Errore critico]: {e}")

# --- CALCOLO STATISTICHE ---
if times:
    avg_time = sum(times) / len(times)
    print("-" * 30)
    print(f"Totale esecuzioni valide: {len(times)}")
    print(f"Tempo Medio: {avg_time:.5f} s")
    print("-" * 30)
else:
    print("\nNessun tempo registrato. Verifica che il comando funzioni.")