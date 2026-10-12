import numpy as np
from scipy.linalg import eigh_tridiagonal
import matplotlib.pyplot as plt

hbar = 1.054571817e-34       # Reduced Planck constant (J s)
mass = 9.1093837139e-31      # Electron mass (kg)
L = 1.0e-9                   # Box length (m)
eV = 1.602176634e-19         # One electronvolt (J)

#Ground state energy for a particle in a box

E_exact = np.pi**2 * hbar**2 / (2 * mass * L**2)

E_exact_eV = E_exact / eV

#Grid sizes for convergence analysis
grid_sizes = [20, 50, 100, 200, 500, 1000]

# Store the calculated results
numerical_energies = []
percentage_errors = []

#Convergence analysis loop

print("=" * 70)
print("PARTICLE IN A BOX: CONVERGENCE ANALYSIS")
print("=" * 70)

print(f"Box length: {L:.3e} m")
print(f"Analytical ground-state energy: {E_exact_eV:.8f} eV")

print("\nResults")
print("-" * 70)

print(
    f"{'Intervals':<15}"
    f"{'Energy (eV)':<20}"
    f"{'Absolute Error (%)':<22}"
)

for N in grid_sizes:

    # Grid spacing and number of interior points
    h = L / N
    M = N - 1

    # Construct the Hamiltonian
    C = hbar**2 / (2 * mass * h**2)

    diagonal = np.full(M, 2 * C)
    off_diagonal = np.full(M - 1, -C)

    # Calculate only the smallest eigenvalue
    energy = eigh_tridiagonal(
        diagonal,
        off_diagonal,
        eigvals_only=True,
        select="i",
        select_range=(0, 0)
    )

    # Ground-state energy in joules
    E_numerical = energy[0]

    # Convert to electronvolts
    E_numerical_eV = E_numerical / eV

    # Calculate percentage error
    error = (
        abs(E_numerical - E_exact) / E_exact
    ) * 100

    # Save the results
    numerical_energies.append(E_numerical_eV)
    percentage_errors.append(error)

    # Print the results
    print(
        f"{N:<15}"
        f"{E_numerical_eV:<20.8f}"
        f"{error:<22.8f}"
    )

print("-" * 70)

print("\nConvergence analysis completed successfully.")

plt.figure(figsize=(8, 5))

plt.plot(
    grid_sizes,
    numerical_energies,
    marker="o",
    label="Numerical ground-state energy"
)

plt.axhline(
    y=E_exact_eV,
    linestyle="--",
    label="Analytical ground-state energy"
)

plt.xlabel("Number of intervals (N)")
plt.ylabel("Ground-state energy (eV)")
plt.title("Ground-State Energy Convergence")
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig("results/energy_convergence.png", dpi=300)
plt.show()

#Percentage error convergence plot

plt.figure(figsize=(8, 5))

plt.loglog(
    grid_sizes,
    percentage_errors,
    marker="o"
)

plt.xlabel("Number of intervals (N)")
plt.ylabel("Absolute percentage error (%)")
plt.title("Numerical Error vs Grid Resolution")
plt.grid(True, which="both")

plt.tight_layout()
plt.savefig("results/error_convergence.png", dpi=300)
plt.show()