# Particle in a Box: Computational Quantum Mechanics

## Overview

This project compares the analytical solution of the one-dimensional infinite potential well with a numerical solution obtained using a finite-difference Hamiltonian. Python, NumPy, SciPy, and Matplotlib are used to calculate energy levels, inspect wavefunctions and probability densities, and study grid convergence.

The model is an electron confined to a box of length \(L=1.0\,\mathrm{nm}\).

## Theory

The potential is

\[
V(x)=
\begin{cases}
0,&0<x<L,\\
\infty,&\text{otherwise}.
\end{cases}
\]

The boundary conditions are \(\psi(0)=\psi(L)=0\). The analytical energies and normalised wavefunctions are

\[
E_n=\frac{n^2\pi^2\hbar^2}{2mL^2},\qquad
\psi_n(x)=\sqrt{\frac{2}{L}}\sin\left(\frac{n\pi x}{L}\right).
\]

The probability density is \(|\psi_n(x)|^2\).

## Numerical method

For a uniform grid with spacing \(h=L/N\), the central finite-difference approximation is

\[
\frac{d^2\psi}{dx^2}\bigg|_{x_i}\approx
\frac{\psi_{i+1}-2\psi_i+\psi_{i-1}}{h^2}.
\]

This produces a symmetric tridiagonal Hamiltonian. Its eigenvalues approximate the energy levels and its eigenvectors represent wavefunction values on the grid. With \(N\) intervals and both boundaries included, there are \(N+1\) grid points and \(N-1\) unknown interior values.

For \(N=200\), the Hamiltonian is \(199\times199\).

## Requirements

- Python 3
- NumPy
- SciPy
- Matplotlib

Activate your virtual environment and install dependencies:

```bash
python -m pip install -r requirements.txt
```

## Run the scripts

Run these commands from the `01_Particle_in_a_Box` project directory:

```bash
python src/analytical_solution.py
python src/numerical_solution.py
python src/convergence_analysis.py
python src/wavefunction_analysis.py
python src/wavefunction_validation.py
python src/final_analysis.py
```

Run from the project root so relative paths to `results/` work correctly.

## Main result

For an electron in a \(1\,\mathrm{nm}\) box, the analytical ground-state energy is approximately \(0.376030\,\mathrm{eV}\). With 200 intervals, the numerical ground-state energy is approximately \(0.376022\,\mathrm{eV}\). Numerical energies approach the analytical values as the grid is refined.

The \(n\)-th state has \(n-1\) internal nodes. The probability density is non-negative and integrates to one for a normalised state.

## Project structure

```text
01_Particle_in_a_Box/
├── src/
│   ├── analytical_solution.py
│   ├── numerical_solution.py
│   ├── convergence_analysis.py
│   ├── wavefunction_analysis.py
│   ├── wavefunction_validation.py
│   └── final_analysis.py
├── notebooks/
├── results/
├── report/
├── requirements.txt
└── README.md
```

## Limitations

- The model assumes an ideal infinite potential well in one dimension.
- Finite differences introduce discretisation error in the energies.
- The agreement between numerical and analytical wavefunctions at grid points is unusually strong for this uniform-grid infinite-well problem; it should not be generalised to arbitrary potentials.
- This is a foundational model, not a molecular electronic-structure calculation.

## Possible extensions

A finite potential well, quantum harmonic oscillator, and hydrogen atom are natural extensions.