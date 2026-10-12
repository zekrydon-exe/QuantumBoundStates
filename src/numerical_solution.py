import numpy as np
from scipy.linalg import eigh_tridiagonal

hbar = 1.054571817e-34       # Reduced Planck constant (J s)
mass = 9.1093837139e-31      # Electron mass (kg)
L = 1.0e-9                   # Box length (m)
eV = 1.602176634e-19         # One electronvolt (J)

N = 200                      # Number of intervals
h = L / N                    # Grid spacing (m)
M = N - 1                    # Number of interior grid points

C = hbar**2 / (2 * mass * h**2)

diagonal = np.full(M, 2 * C)
off_diagonal = np.full(M - 1, -C)

#Eigenvalue problem solution using the tridiagonal matrix algorithm

energy, wavefunctions = eigh_tridiagonal(
    diagonal,
    off_diagonal
)

# Convert numerical energies from joules to electronvolts
energy_eV = energy / eV

n = np.arange(1, 6)

analytical_energy = (
    n**2 * np.pi**2 * hbar**2
) / (2 * mass * L**2)

analytical_energy_eV = analytical_energy / eV

absolute_error = np.abs(
    energy_eV[:5] - analytical_energy_eV
)

percentage_error = (
    absolute_error / analytical_energy_eV
) * 100

print("=" * 78)
print("PARTICLE IN A BOX: NUMERICAL SOLUTION")
print("=" * 78)

print(f"Box length             : {L:.3e} m")
print(f"Number of intervals    : {N}")
print(f"Number of interior points: {M}")

print("\nEnergy Comparison")
print("-" * 78)

print(
    f"{'n':<5}"
    f"{'Analytical (eV)':<22}"
    f"{'Numerical (eV)':<22}"
    f"{'Error (%)':<15}"
)

for i in range(5):
    print(
        f"{i + 1:<5}"
        f"{analytical_energy_eV[i]:<22.8f}"
        f"{energy_eV[i]:<22.8f}"
        f"{percentage_error[i]:<15.6f}"
    )

print("-" * 78)

print(f"\nAnalytical ground-state energy : {analytical_energy_eV[0]:.8f} eV")
print(f"Numerical ground-state energy  : {energy_eV[0]:.8f} eV")
print(f"Ground-state percentage error  : {percentage_error[0]:.6f}%")

print("\nCalculation completed successfully.")