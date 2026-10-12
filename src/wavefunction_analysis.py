import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal
from pathlib import Path

hbar = 1.054571817e-34       # Reduced Planck constant (J s)
mass = 9.1093837139e-31      # Electron mass (kg)
L = 1.0e-9                   # Box length (m)
eV = 1.602176634e-19         # One electronvolt (J)

N = 200
h = L / N
M = N - 1

# Interior grid points
x = np.arange(1, N) * h

# Including the boundaries for plotting
x_full = np.linspace(0, L, N + 1)

# Dimensionless position: x/L
x_dimensionless = x / L
x_full_dimensionless = x_full / L

#Hamiltonian matrix construction

C = hbar**2 / (2 * mass * h**2)

diagonal = np.full(M, 2 * C)
off_diagonal = np.full(M - 1, -C)

# Calculate eigenvalues and eigenvectors
energy, eigenvectors = eigh_tridiagonal(
    diagonal,
    off_diagonal
)

energy_eV = energy / eV

number_of_states = 4

psi_numerical = eigenvectors[:, :number_of_states] / np.sqrt(h)

#Calculation of Analytical Wavefunctions

psi_analytical = np.zeros(
    (N + 1, number_of_states)
)

for i in range(number_of_states):

    n = i + 1

    psi_analytical[:, i] = (
        np.sqrt(2 / L)
        * np.sin(n * np.pi * x_full / L)
    )

# Add boundary values to numerical wavefunctions
psi_numerical_full = np.zeros(
    (N + 1, number_of_states)
)

psi_numerical_full[1:N, :] = psi_numerical

for i in range(number_of_states):

    overlap = np.sum(
        psi_numerical_full[:, i]
        * psi_analytical[:, i]
    ) * h

    if overlap < 0:
        psi_numerical_full[:, i] *= -1


results_directory = Path("results")
results_directory.mkdir(parents=True, exist_ok=True)

#Plotting of Wavefunctions

plt.figure(figsize=(10, 6))

for i in range(number_of_states):

    n = i + 1

    # Scale the wavefunction by sqrt(L) to make it dimensionless
    numerical_scaled = (
        np.sqrt(L) * psi_numerical_full[:, i]
    )

    analytical_scaled = (
        np.sqrt(L) * psi_analytical[:, i]
    )

    plt.plot(
        x_full_dimensionless,
        numerical_scaled,
        label=f"Numerical n={n}"
    )

    plt.plot(
        x_full_dimensionless,
        analytical_scaled,
        linestyle="--",
        label=f"Analytical n={n}"
    )

plt.xlabel("Position, x/L")
plt.ylabel(r"Dimensionless wavefunction, $\sqrt{L}\,\psi(x)$")
plt.title("Particle in a Box: Numerical vs Analytical Wavefunctions")
plt.legend(ncol=2)
plt.grid(True)
plt.tight_layout()

plt.savefig(
    results_directory / "wavefunctions.png",
    dpi=300
)

plt.show()

#Plotting of Probability Densities

plt.figure(figsize=(10, 6))

for i in range(number_of_states):

    n = i + 1

    # L * |psi|^2 is a dimensionless probability density
    numerical_density = (
        L * np.abs(psi_numerical_full[:, i])**2
    )

    analytical_density = (
        L * np.abs(psi_analytical[:, i])**2
    )

    plt.plot(
        x_full_dimensionless,
        numerical_density,
        label=f"Numerical n={n}"
    )

    plt.plot(
        x_full_dimensionless,
        analytical_density,
        linestyle="--",
        label=f"Analytical n={n}"
    )

plt.xlabel("Position, x/L")
plt.ylabel(r"Dimensionless probability density, $L|\psi(x)|^2$")
plt.title("Particle in a Box: Probability Densities")
plt.legend(ncol=2)
plt.grid(True)
plt.tight_layout()

plt.savefig(
    results_directory / "probability_densities.png",
    dpi=300
)

plt.show()

#Display Energy Levels

print("=" * 65)
print("PARTICLE IN A BOX: WAVEFUNCTION ANALYSIS")
print("=" * 65)

print(f"Box length: {L:.3e} m")
print(f"Grid intervals: {N}")

print("\nFirst four numerical energy levels:")

for i in range(number_of_states):
    print(
        f"n = {i + 1}: "
        f"E = {energy_eV[i]:.8f} eV"
    )

print("\nPlots saved:")
print("results/wavefunctions.png")
print("results/probability_densities.png")

print("\nWavefunction analysis completed successfully.")