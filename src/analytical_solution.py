import numpy as np

from quantum_well import analytical_energies, ELECTRONVOLT

L = 1.0e-9  # Box length (m)
n = np.arange(1, 6)

energy = analytical_energies(L, n)
energy_eV = energy / ELECTRONVOLT

print("Particle in a 1D Infinite Potential Well")
print(f"Box length: {L:.2e} m")
print()

for quantum_number, E in zip(n, energy_eV):
    print(f"n = {quantum_number}: E = {E:.6f} eV")