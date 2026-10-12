"""Shared calculations for a particle in a 1D infinite potential well."""

import numpy as np
from scipy.linalg import eigh_tridiagonal

HBAR = 1.054571817e-34       # Reduced Planck constant (J s)
ELECTRON_MASS = 9.1093837139e-31  # Electron mass (kg)
ELECTRONVOLT = 1.602176634e-19    # One electronvolt (J)

#Analytical solution for a particle in a 1D infinite potential well
def analytical_energies(L, quantum_numbers):
    """Calculate analytical energy eigenvalues in joules.

    Parameters
    ----------
    L : float
        Length of the box in metres.
    quantum_numbers : array-like
        Positive integer quantum numbers.

    Returns
    -------
    numpy.ndarray
        Analytical energy eigenvalues in joules.
    """
    n = np.asarray(quantum_numbers, dtype=int)

    if L <= 0:
        raise ValueError("Box length L must be positive.")

    if np.any(n < 1):
        raise ValueError("Quantum numbers must be positive integers.")

    return (
        n**2 * np.pi**2 * HBAR**2
        / (2.0 * ELECTRON_MASS * L**2)
    )

#Numerical solution for a particle in a 1D infinite potential well
def numerical_solution(L, N, number_of_states=5):
    """Solve the finite-difference Hamiltonian.

    Parameters
    ----------
    L : float
        Length of the box in metres.
    N : int
        Number of grid intervals.
    number_of_states : int
        Number of lowest-energy states to calculate.

    Returns
    -------
    x : numpy.ndarray
        Full position grid, including both boundaries.
    energy : numpy.ndarray
        Numerical energies in joules.
    wavefunctions : numpy.ndarray
        Continuum-scaled numerical wavefunctions at interior grid points.
        Shape: (N - 1, number_of_states).
    """
    if L <= 0:
        raise ValueError("Box length L must be positive.")

    if not isinstance(N, (int, np.integer)) or N < 2:
        raise ValueError("N must be an integer greater than or equal to 2.")

    if (
        not isinstance(number_of_states, (int, np.integer))
        or not 1 <= number_of_states <= N - 1
    ):
        raise ValueError(
            "number_of_states must be between 1 and N - 1."
        )

    h = L / N
    interior_points = N - 1

    coefficient = HBAR**2 / (
        2.0 * ELECTRON_MASS * h**2
    )

    diagonal = np.full(interior_points, 2.0 * coefficient)
    off_diagonal = np.full(interior_points - 1, -coefficient)

    energy, eigenvectors = eigh_tridiagonal(
        diagonal,
        off_diagonal,
        select="i",
        select_range=(0, number_of_states - 1),
    )

    # Convert discrete normalisation to approximate continuum normalisation.
    wavefunctions = eigenvectors / np.sqrt(h)

    x = np.linspace(0.0, L, N + 1)

    return x, energy, wavefunctions