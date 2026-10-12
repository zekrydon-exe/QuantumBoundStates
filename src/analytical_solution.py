import numpy as np

# Physical constants (SI units)
hbar = 1.054571817e-34  # Reduced Planck constant, J s
mass = 9.1093837139e-31  # Electron mass, kg
L = 1.0e-9  # Box length, m

# Quantum numbers n = 1, 2, 3, 4, 5
n = np.arange(1, 6)

# Calculate energy levels in joules
energy = (n**2 * np.pi**2 * hbar**2) / (2 * mass * L**2)

# Convert joules to electronvolts
energy_eV = energy / 1.602176634e-19

# Display results
print("Particle in a 1D Infinite Potential Well")
print(f"Box length: {L:.2e} m")
print()

for quantum_number, E in zip(n, energy_eV):
    print(f"n = {quantum_number}: E = {E:.6f} eV")