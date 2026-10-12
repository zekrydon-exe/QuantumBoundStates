import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.integrate import trapezoid

hbar = 1.054571817e-34
mass = 9.1093837139e-31
L = 1.0e-9
eV = 1.602176634e-19

#Numberical grid parameters

N = 200
h = L / N
M = N - 1

x = np.linspace(0, L, N + 1)

# Interior points
x_interior = x[1:-1]

C = hbar**2 / (2 * mass * h**2)

diagonal = np.full(M, 2 * C)
off_diagonal = np.full(M - 1, -C)

energy, eigenvectors = eigh_tridiagonal(
    diagonal,
    off_diagonal,
    select="i",
    select_range=(0, 4)
)

energy_eV = energy / eV

psi_numerical = eigenvectors / np.sqrt(h)

# Add boundary values
psi_full = np.zeros((N + 1, 5))
psi_full[1:-1, :] = psi_numerical

#Analytical wavefunctions for comparison

psi_exact = np.zeros((N + 1, 5))

for i in range(5):

    n = i + 1

    psi_exact[:, i] = (
        np.sqrt(2 / L)
        * np.sin(n * np.pi * x / L)
    )

#Normalization of analytical wavefunctions

print("=" * 85)
print("WAVEFUNCTION VALIDATION")
print("=" * 85)

print(
    f"{'State':<10}"
    f"{'Energy (eV)':<18}"
    f"{'Norm':<18}"
    f"{'Relative L2 error':<22}"
)

for i in range(5):

    # Numerical wavefunction
    psi_num = psi_full[:, i].copy()

    # Analytical wavefunction
    psi_ref = psi_exact[:, i]

    # Align arbitrary eigenvector sign
    overlap = trapezoid(
        psi_num * psi_ref,
        x
    )

    if overlap < 0:
        psi_num *= -1

    # Continuous normalisation integral
    norm = trapezoid(
        np.abs(psi_num)**2,
        x
    )

    # Relative L2 error
    difference = psi_num - psi_ref

    numerator = np.sqrt(
        trapezoid(np.abs(difference)**2, x)
    )

    denominator = np.sqrt(
        trapezoid(np.abs(psi_ref)**2, x)
    )

    relative_error = numerator / denominator

    print(
        f"{i + 1:<10}"
        f"{energy_eV[i]:<18.8f}"
        f"{norm:<18.10f}"
        f"{relative_error:<22.8e}"
    )

print("=" * 85)
print("Validation completed.")