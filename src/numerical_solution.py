"""Calculate and compare analytical and numerical energy levels."""

import numpy as np

from quantum_well import (
    analytical_energies,
    numerical_solution,
    ELECTRONVOLT,
)

L = 1.0e-9       # Box length (m)
N = 200          # Number of grid intervals
number_of_states = 5

n = np.arange(1, number_of_states + 1)

analytical_energy = analytical_energies(L, n)
analytical_energy_eV = analytical_energy / ELECTRONVOLT

x, energy, wavefunctions = numerical_solution(
    L,
    N,
    number_of_states,
)

energy_eV = energy / ELECTRONVOLT

#Energy error calculations
absolute_error = np.abs(
    energy_eV - analytical_energy_eV
)

percentage_error = (
    absolute_error / analytical_energy_eV
) * 100.0

print("=" * 78)
print("PARTICLE IN A BOX: NUMERICAL SOLUTION")
print("=" * 78)

print(f"Box length               : {L:.3e} m")
print(f"Number of intervals      : {N}")
print(f"Number of interior points: {N - 1}")

print("\nEnergy Comparison")
print("-" * 78)

print(
    f"{'n':<5}"
    f"{'Analytical (eV)':<22}"
    f"{'Numerical (eV)':<22}"
    f"{'Error (%)':<15}"
)

for i in range(number_of_states):
    print(
        f"{n[i]:<5}"
        f"{analytical_energy_eV[i]:<22.8f}"
        f"{energy_eV[i]:<22.8f}"
        f"{percentage_error[i]:<15.6f}"
    )

print("-" * 78)

print(
    f"\nAnalytical ground-state energy : "
    f"{analytical_energy_eV[0]:.8f} eV"
)

print(
    f"Numerical ground-state energy  : "
    f"{energy_eV[0]:.8f} eV"
)

print(
    f"Ground-state percentage error  : "
    f"{percentage_error[0]:.6f}%"
)

print("\nCalculation completed successfully.")