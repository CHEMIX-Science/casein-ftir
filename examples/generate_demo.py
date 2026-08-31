"""Generate synthetic demonstrations, never experimental measurements."""
from pathlib import Path
from casein_ftir.simulate import simulate_casein, simulate_galalithe
from casein_ftir.io_module import save_spectrum_csv

def main():
    folder = Path(__file__).resolve().parent
    casein = simulate_casein(noise_level=0.0)
    galalithe = simulate_galalithe(crosslink_degree=0.7, residual_formol=0.02, noise_level=0.0)
    save_spectrum_csv(casein, folder / "synthetic_casein.csv")
    save_spectrum_csv(galalithe, folder / "synthetic_galalithe.csv")
    print("Generated two synthetic demo CSVs; not experimental data.")

if __name__ == "__main__":
    main()
