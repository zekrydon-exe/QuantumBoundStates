# Particle in a Box: Computational Quantum Mechanics

## Overview

This project investigates the one-dimensional infinite potential well using analytical quantum mechanics and numerical computation.

The analytical solution is compared with a numerical solution obtained by discretising the time-independent Schrödinger equation using the finite-difference method.

**Key objectives:**
- Calculate quantised energy levels analytically and numerically.
- Construct and diagonalise the Hamiltonian matrix.
- Visualise wavefunctions and probability densities.
- Investigate the effect of grid resolution on numerical accuracy.
- Validate numerical results against analytical solutions.

**System:** An electron confined to a one-dimensional box of length \(L = 1.0\,\mathrm{nm}\).

## 1. Theoretical Background

### 1.1 Infinite Potential Well

The potential energy is defined as

$$
V(x) =
\begin{cases}
0, & 0 < x < L,\\
\infty, & \text{otherwise}.
\end{cases}
$$

Inside the box, the time-independent Schrödinger equation becomes

$$
-\frac{\hbar^2}{2m}\frac{d^2\psi(x)}{dx^2}
= E\psi(x).
$$

The infinite potential barriers impose the boundary conditions

$$
\psi(0) = \psi(L) = 0.
$$

### 1.2 Analytical Energy Levels

Applying the boundary conditions gives the allowed energy eigenvalues

$$
E_n = \frac{n^2\pi^2\hbar^2}{2mL^2},
\qquad n = 1,2,3,\ldots
$$

The ground-state energy is non-zero, and the energy levels scale as

$$
E_n \propto \frac{n^2}{L^2}.
$$

Consequently, doubling the box length reduces every energy level by a factor of four.

### 1.3 Analytical Wavefunctions

The normalised stationary-state wavefunctions are

$$
\psi_n(x) =
\sqrt{\frac{2}{L}}
\sin\left(\frac{n\pi x}{L}\right).
$$

The probability density is

$$
P_n(x) = |\psi_n(x)|^2.
$$

Normalisation requires

$$
\int_0^L |\psi_n(x)|^2\,dx = 1.
$$

The \(n\)-th state contains \(n-1\) internal nodes.

## 2. Numerical Method

### 2.1 Spatial Discretisation

The box is divided into \(N\) equal intervals with grid spacing

$$
h = \frac{L}{N}.
$$

The grid points are

$$
x_i = ih,
\qquad i=0,1,\ldots,N.
$$

Both boundaries are included, giving \(N+1\) grid points and \(N-1\) unknown interior wavefunction values.

### 2.2 Finite-Difference Approximation

The second derivative is approximated using the central finite-difference formula

$$
\frac{d^2\psi}{dx^2}\bigg|_{x_i}
\approx
\frac{\psi_{i+1}-2\psi_i+\psi_{i-1}}{h^2}.
$$

Substituting this approximation into the Schrödinger equation produces a matrix eigenvalue problem

$$
\mathbf{H}\boldsymbol{\psi}_n
= E_n\boldsymbol{\psi}_n.
$$

The resulting Hamiltonian is a symmetric tridiagonal matrix.

Its diagonal elements are

$$
H_{ii} = \frac{\hbar^2}{mh^2},
$$

and its neighbouring off-diagonal elements are

$$
H_{i,i+1} = H_{i+1,i}
= -\frac{\hbar^2}{2mh^2}.
$$

The eigenvalues approximate the allowed energies, while the eigenvectors represent the discretised wavefunctions.

### 2.3 Numerical Implementation

The project uses Python and the following libraries:

- **NumPy:** Numerical arrays and mathematical operations.
- **SciPy:** Eigenvalue calculations using a specialised tridiagonal solver.
- **Matplotlib:** Visualisation of energies, wavefunctions, probability densities and convergence.

The numerical energies are calculated independently of the analytical energy formula. The analytical results are subsequently used to assess numerical accuracy.

## 3. Results and Validation

### 3.1 Energy Levels

For an electron confined to a box of length \(1\,\mathrm{nm}\), the analytical ground-state energy is approximately

$$
E_1 \approx 0.376030\,\mathrm{eV}.
$$

With 200 grid intervals, the numerical ground-state energy is approximately

$$
E_1^{\mathrm{num}} \approx 0.376022\,\mathrm{eV}.
$$

The small difference arises from the finite-difference approximation.

See `results/energy_comparison.png` for the energy-level comparison.

### 3.2 Grid Convergence

The numerical energy is calculated for several grid resolutions.

For sufficiently fine grids, the leading finite-difference energy error scales approximately as

$$
\epsilon \propto h^2
\propto \frac{1}{N^2}.
$$

Increasing the number of intervals generally improves the energy estimate, although it also increases computational cost.

See `results/energy_convergence.png` for the convergence analysis.

### 3.3 Wavefunctions and Probability Densities

The numerical wavefunctions are compared with the analytical solutions.

The probability density is non-negative, and the normalisation condition requires its integral over the box to equal one.

The sampled numerical and analytical wavefunctions agree exceptionally closely for this particular uniform-grid infinite-well problem. This is a special mathematical property of the model and should not be generalised to arbitrary potentials.

See the following figures:

- `results/wavefunctions.png`
- `results/probability_densities.png`

### 3.4 Error Analysis

The relative percentage error in an energy eigenvalue is calculated as

$$
\epsilon_n(\%) =
\frac{
|E_n^{\mathrm{num}}-E_n^{\mathrm{exact}}|
}{
|E_n^{\mathrm{exact}}|
}
\times 100.
$$

The numerical errors and convergence data are saved as CSV files for further analysis.

## 4. Project Structure

```text
Particle-in-a-Box/
├── notebooks/
│   └── exploration.ipynb
├── report/
│   └── technical_report.md
├── results/
│   ├── convergence_data.csv
│   ├── energy_comparison.csv
│   ├── energy_comparison.png
│   ├── energy_convergence.png
│   ├── energy_errors.png
│   ├── error_convergence.png
│   ├── probability_densities.png
│   └── wavefunctions.png
├── src/
│   ├── analytical_solution.py
│   ├── final_analysis.py
│   ├── numerical_solution.py
│   ├── quantum_well.py
│   └── wavefunction_validation.py
├── .gitignore
├── README.md
└── requirements.txt
```

## 5. Code Structure
The project is organised into reusable calculation modules, standalone scripts, and analysis tools. The shared module, quantum_well.py, contains the core analytical and numerical calculations used throughout the project.

1. src/quantum_well.py - Defines physical constants and reusable functions for calculating analytical energy eigenvalues and solving the finite-difference Hamiltonian.
2. src/analytical_solution.py - Calculates and displays the analytical energy eigenvalues for a particle in a one-dimensional infinite potential well.
3. src/numerical_solution.py - Calculates numerical energy eigenvalues using the finite-difference method and compares them with the analytical results.
4. src/final_analysis.py - Performs the final numerical analysis, generates energy comparisons, plots wavefunctions and probability densities, and investigates convergence with grid refinement. Saves plots and numerical data in the results/ directory.
5. src/wavefunction_validation.py - Checks wavefunction normalisation and calculates relative L2 errors to assess agreement between numerical and analytical wavefunctions.

## 6. Installation

Python 3 is required.

Install the dependencies in your activated virtual environment:

```bash
python -m pip install -r requirements.txt
```

## 7. Running the Project

Run the commands from the project root directory:

```bash
python src/analytical_solution.py
python src/numerical_solution.py
python src/wavefunction_validation.py
python src/final_analysis.py
```

The scripts generate numerical results, plots and CSV files in the `results/` directory.

## 8. Limitations

- The model assumes an ideal one-dimensional infinite potential well.
- Finite-difference discretisation introduces numerical error.
- The sampled wavefunction comparison is unusually favourable for this particular system.
- The model does not represent a realistic molecular electronic structure calculation.

## 9. Future Extensions

Potential extensions include:

- Finite potential wells.
- The quantum harmonic oscillator.
- The hydrogen atom and hydrogenic orbitals.
- More general one-dimensional potentials.
- Additional numerical methods for solving the Schrödinger equation.

## Author

Computational quantum mechanics project exploring the relationship between analytical quantum theory and numerical methods.