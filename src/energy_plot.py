import matplotlib.pyplot as plt
import csv
import math

def load_energy_data(filename):
    times = []
    energies = []
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            try:
                t = float(row[0])
                e = float(row[1])
                if not (math.isfinite(e)):
                    continue
                times.append(t)
                energies.append(e)
            except ValueError:
                continue
    return times, energies


times_rk4, energy_rk4 = load_energy_data("energy_rk4.csv")
times_euler, energy_euler = load_energy_data("energy_euler.csv")
times_sym, energy_sym = load_energy_data("energy_sym.csv")
times_verlet, energy_verlet = load_energy_data("energy_ver.csv")

e0_rk4 = energy_rk4[0]
e0_euler = energy_euler[0]
e0_sym = energy_sym[0]
e0_verlet = energy_verlet[0]

error_rk4 = [e - e0_rk4 for e in energy_rk4]
error_euler = [e - e0_euler for e in energy_euler]
error_sym = [e - e0_sym for e in energy_sym]
error_verlet = [e - e0_verlet for e in energy_verlet]

plt.figure(figsize=(10,6))
plt.plot(times_rk4, error_rk4, label="RK4", linewidth=2)
plt.plot(times_euler, error_euler, label="Euler", linewidth=2)
plt.plot(times_sym, error_sym, label="Symplectic", linewidth=2)
plt.plot(times_verlet, error_verlet, label="Verlet", linewidth=2)

plt.xlabel("Time")
plt.ylabel("Energy Error (E(t) - E₀)")
plt.title("Energy Error Over Time for Different Integrators")
plt.legend()
plt.grid(True)
plt.savefig("energy_error_comparison.png", dpi=300)
plt.show()
