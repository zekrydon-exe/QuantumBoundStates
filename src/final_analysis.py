import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from quantum_well import (
    HBAR,
    ELECTRON_MASS,
    ELECTRONVOLT,
    analytical_energies,
    numerical_solution,
)


L = 1.0e-9
N = 200
number_of_states = 5

PROJECT_ROOT = Path(__file__).resolve().parent.parent
results_dir = PROJECT_ROOT / "results"
results_dir.mkdir(parents=True, exist_ok=True)

n = np.arange(1, number_of_states + 1)
x, energy, psi_num = numerical_solution(L, N, number_of_states)
E_exact = analytical_energies(L, n)

E_exact_eV = E_exact / ELECTRONVOLT
E_num_eV = energy / ELECTRONVOLT

M = N - 1
psi_num_full = np.zeros((N + 1, number_of_states))
psi_num_full[1:-1, :] = psi_num

psi_exact = np.zeros((N + 1, number_of_states))
for j, quantum_number in enumerate(n):
    psi_exact[:, j] = (
        np.sqrt(2.0 / L)
        * np.sin(quantum_number * np.pi * x / L)
    )

for j in range(number_of_states):
    overlap = np.trapezoid(psi_num_full[:, j] * psi_exact[:, j], x)
    if overlap < 0.0:
        psi_num_full[:, j] *= -1.0

absolute_energy_error = np.abs(E_num_eV - E_exact_eV)
percentage_energy_error = absolute_energy_error / E_exact_eV * 100.0

print("=" * 85)
print("PARTICLE IN A BOX: FINAL NUMERICAL ANALYSIS")
print("=" * 85)
print(f"Box length: {L:.4e} m")
print(f"Grid intervals: {N}")
print(f"Interior grid points: {M}")

print("\nEnergy comparison")
print("-" * 85)
print(
    f"{'n':<6}"
    f"{'Analytical (eV)':<22}"
    f"{'Numerical (eV)':<22}"
    f"{'Error (%)':<18}"
)

for j in range(number_of_states):
    print(
        f"{n[j]:<6}"
        f"{E_exact_eV[j]:<22.10f}"
        f"{E_num_eV[j]:<22.10f}"
        f"{percentage_energy_error[j]:<18.8e}"
    )

print("-" * 85)

energy_data = np.column_stack(
    (n, E_exact_eV, E_num_eV, absolute_energy_error, percentage_energy_error)
)
np.savetxt(
    results_dir / "energy_comparison.csv",
    energy_data,
    delimiter=",",
    header=(
        "quantum_number,analytical_energy_eV,numerical_energy_eV,"
        "absolute_error_eV,percentage_error"
    ),
    comments="",
)

plt.figure(figsize=(9, 5))
plt.plot(n, E_exact_eV, marker="o", linestyle="-", label="Analytical")
plt.plot(n, E_num_eV, marker="x", linestyle="--", label="Numerical")
plt.xlabel("Quantum number, n")
plt.ylabel("Energy (eV)")
plt.title("Analytical vs Numerical Energy Levels")
plt.xticks(n)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(results_dir / "energy_comparison.png", dpi=300)
plt.close()

plt.figure(figsize=(9, 5))
plt.semilogy(n, percentage_energy_error, marker="o")
plt.xlabel("Quantum number, n")
plt.ylabel("Absolute percentage error (%)")
plt.title("Energy Error for the First Five States")
plt.xticks(n)
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig(results_dir / "energy_errors.png", dpi=300)
plt.close()

plt.figure(figsize=(10, 6))
for j in range(number_of_states):
    numerical_scaled = np.sqrt(L) * psi_num_full[:, j]
    analytical_scaled = np.sqrt(L) * psi_exact[:, j]
    plt.plot(x / L, numerical_scaled, label=f"Numerical n={j + 1}")
    plt.plot(
        x / L,
        analytical_scaled,
        linestyle="--",
        label=f"Analytical n={j + 1}",
    )

plt.xlabel("Dimensionless position, x/L")
plt.ylabel(r"$\sqrt{L}\,\psi(x)$")
plt.title("Numerical and Analytical Wavefunctions")
plt.legend(ncol=2)
plt.grid(True)
plt.tight_layout()
plt.savefig(results_dir / "wavefunctions.png", dpi=300)
plt.close()

plt.figure(figsize=(10, 6))
for j in range(number_of_states):
    numerical_density = L * np.abs(psi_num_full[:, j]) ** 2
    analytical_density = L * np.abs(psi_exact[:, j]) ** 2
    plt.plot(x / L, numerical_density, label=f"Numerical n={j + 1}")
    plt.plot(
        x / L,
        analytical_density,
        linestyle="--",
        label=f"Analytical n={j + 1}",
    )

plt.xlabel("Dimensionless position, x/L")
plt.ylabel(r"$L|\psi(x)|^2$")
plt.title("Numerical and Analytical Probability Densities")
plt.legend(ncol=2)
plt.grid(True)
plt.tight_layout()
plt.savefig(results_dir / "probability_densities.png", dpi=300)
plt.close()

grid_sizes = np.array([20, 50, 100, 200, 500, 1000], dtype=int)
convergence_energies = []
convergence_errors = []

for grid_N in grid_sizes:
    grid_h = L / grid_N
    grid_M = grid_N - 1
    grid_C = HBAR**2 / (2.0 * ELECTRON_MASS * grid_h**2)

    grid_diagonal = np.full(grid_M, 2.0 * grid_C)
    grid_off_diagonal = np.full(grid_M - 1, -grid_C)

    ground_energy = numerical_solution(
        L, grid_N, number_of_states=1
    )[1][0]

    convergence_energies.append(ground_energy / ELECTRONVOLT)
    convergence_errors.append(
        abs(ground_energy - E_exact[0]) / E_exact[0] * 100.0
    )

convergence_energies = np.asarray(convergence_energies)
convergence_errors = np.asarray(convergence_errors)

convergence_data = np.column_stack(
    (grid_sizes, convergence_energies, convergence_errors)
)
np.savetxt(
    results_dir / "convergence_data.csv",
    convergence_data,
    delimiter=",",
    header="intervals,ground_energy_eV,percentage_error",
    comments="",
)

print("\nGround-state convergence")
print("-" * 65)
print(f"{'Intervals':<15}{'Energy (eV)':<22}{'Error (%)':<18}")
for i, grid_N in enumerate(grid_sizes):
    print(
        f"{grid_N:<15}"
        f"{convergence_energies[i]:<22.10f}"
        f"{convergence_errors[i]:<18.8e}"
    )
print("-" * 65)

plt.figure(figsize=(9, 5))
plt.plot(
    grid_sizes,
    convergence_energies,
    marker="o",
    label="Numerical ground state",
)
plt.axhline(
    E_exact_eV[0],
    linestyle="--",
    label="Analytical ground state",
)
plt.xlabel("Number of grid intervals")
plt.ylabel("Ground-state energy (eV)")
plt.title("Ground-State Energy Convergence")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(results_dir / "energy_convergence.png", dpi=300)
plt.close()

plt.figure(figsize=(9, 5))
plt.loglog(grid_sizes, convergence_errors, marker="o")
plt.xlabel("Number of grid intervals")
plt.ylabel("Absolute percentage error (%)")
plt.title("Numerical Error vs Grid Resolution")
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig(results_dir / "error_convergence.png", dpi=300)
plt.close()

print("\nFiles saved in the results directory:")
for file in sorted(results_dir.iterdir()):
    print(f"  {file.name}")

print("\nFinal analysis completed successfully.")