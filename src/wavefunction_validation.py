import numpy as np
from scipy.integrate import trapezoid

from quantum_well import (
    numerical_solution,
    ELECTRONVOLT,
)


L = 1.0e-9
N = 200
number_of_states = 5

NORM_TOLERANCE = 1e-10
L2_ERROR_TOLERANCE = 1e-10

n = np.arange(1, number_of_states + 1)

x, energy, psi_num = numerical_solution(
    L,
    N,
    number_of_states,
)

energy_eV = energy / ELECTRONVOLT

psi_num_full = np.zeros((N + 1, number_of_states))
psi_num_full[1:-1, :] = psi_num

psi_exact = np.zeros((N + 1, number_of_states))

for j, quantum_number in enumerate(n):
    psi_exact[:, j] = (
        np.sqrt(2.0 / L)
        * np.sin(quantum_number * np.pi * x / L)
    )

for j in range(number_of_states):
    overlap = trapezoid(
        psi_num_full[:, j] * psi_exact[:, j],
        x,
    )

    if overlap < 0.0:
        psi_num_full[:, j] *= -1.0

validation_passed = True

print("=" * 85)
print("WAVEFUNCTION VALIDATION")
print("=" * 85)

print(
    f"{'State':<10}"
    f"{'Energy (eV)':<18}"
    f"{'Norm':<18}"
    f"{'Relative L2 error':<22}"
)

for i in range(number_of_states):
    psi = psi_num_full[:, i]
    psi_reference = psi_exact[:, i]

    norm = trapezoid(np.abs(psi) ** 2, x)

    difference = psi - psi_reference

    numerator = np.sqrt(
        trapezoid(np.abs(difference) ** 2, x)
    )

    denominator = np.sqrt(
        trapezoid(np.abs(psi_reference) ** 2, x)
    )

    relative_error = numerator / denominator

    print(
        f"{i + 1:<10}"
        f"{energy_eV[i]:<18.8f}"
        f"{norm:<18.10f}"
        f"{relative_error:<22.8e}"
    )

    if abs(norm - 1.0) > NORM_TOLERANCE:
        validation_passed = False

    if relative_error > L2_ERROR_TOLERANCE:
        validation_passed = False

print("=" * 85)

if validation_passed:
    print("VALIDATION PASSED")
else:
    print("VALIDATION FAILED")
    raise SystemExit(1)