# Technical Report: Numerical Solution of the One-Dimensional Infinite Potential Well

## Abstract

The one-dimensional infinite potential well is a fundamental model in quantum mechanics that demonstrates energy quantisation and the role of boundary conditions in determining allowed quantum states. This project investigates the system using both analytical quantum mechanics and numerical computation.

The time-independent Schrödinger equation is discretised using a central finite-difference approximation, resulting in a symmetric tridiagonal Hamiltonian matrix. The matrix is diagonalised using Python and SciPy to obtain numerical energy eigenvalues and eigenvectors.

The numerical energy levels are compared with their analytical counterparts. Wavefunctions, probability densities, normalisation, and grid convergence are also investigated. For an electron confined to a box of length \(1.0\,\mathrm{nm}\), the analytical ground-state energy is approximately \(0.376030\,\mathrm{eV}\), while the numerical result obtained using 200 grid intervals is approximately \(0.376022\,\mathrm{eV}\).

The results demonstrate how numerical methods can approximate quantum-mechanical systems and how analytical solutions can be used to assess numerical accuracy.

---

## 1. Introduction

Quantum mechanics describes the behaviour of microscopic particles through wavefunctions and operators. Unlike classical mechanics, quantum mechanics often predicts discrete energy values for confined systems.

The one-dimensional infinite potential well, commonly called the particle-in-a-box model, is one of the simplest systems exhibiting energy quantisation.

Although an analytical solution exists for this system, it also provides an excellent benchmark for learning computational quantum mechanics. By comparing numerical results with exact analytical solutions, the accuracy and limitations of numerical methods can be investigated.

### 1.1 Objectives

The primary objectives of this project are:

1. To derive the analytical energy eigenvalues and wavefunctions of a particle confined in a one-dimensional infinite potential well.
2. To discretise the time-independent Schrödinger equation using the finite-difference method.
3. To construct and diagonalise the Hamiltonian matrix numerically.
4. To compare numerical energy eigenvalues with analytical energy levels.
5. To visualise wavefunctions and probability densities.
6. To investigate numerical convergence as the spatial grid is refined.
7. To evaluate the normalisation and accuracy of numerical wavefunctions.

---

## 2. Theoretical Background

### 2.1 The Infinite Potential Well

Consider a particle of mass \(m\) confined to a one-dimensional box of length \(L\).

The potential energy is defined as

$$
V(x)=
\begin{cases}
0, & 0<x<L,\\
\infty, & \text{otherwise}.
\end{cases}
$$

Inside the box, the particle experiences zero potential energy. Outside the box, the potential is infinite, preventing the particle from occupying the external region.

Consequently, the wavefunction must satisfy the boundary conditions

$$
\psi(0)=0
$$

and

$$
\psi(L)=0.
$$

These boundary conditions restrict the allowed wavefunctions and lead to discrete energy levels.

### 2.2 Time-Independent Schrödinger Equation

The one-dimensional time-independent Schrödinger equation is

$$
-\frac{\hbar^2}{2m}
\frac{d^2\psi(x)}{dx^2}
+
V(x)\psi(x)
=
E\psi(x).
$$

Inside the potential well, \(V(x)=0\). Therefore,

$$
-\frac{\hbar^2}{2m}
\frac{d^2\psi(x)}{dx^2}
=
E\psi(x).
$$

Rearranging gives

$$
\frac{d^2\psi(x)}{dx^2}
+
\frac{2mE}{\hbar^2}\psi(x)
=
0.
$$

Defining

$$
k^2=\frac{2mE}{\hbar^2},
$$

the equation becomes

$$
\frac{d^2\psi(x)}{dx^2}+k^2\psi(x)=0.
$$

Its general solution is

$$
\psi(x)=A\sin(kx)+B\cos(kx).
$$

Applying the first boundary condition, \(\psi(0)=0\), gives \(B=0\).

The second boundary condition requires

$$
A\sin(kL)=0.
$$

For a nontrivial wavefunction, \(A\neq 0\), so

$$
\sin(kL)=0.
$$

Therefore,

$$
kL=n\pi,
\qquad n=1,2,3,\ldots
$$

and hence

$$
k_n=\frac{n\pi}{L}.
$$

### 2.3 Analytical Energy Eigenvalues

Using the relationship

$$
E=\frac{\hbar^2k^2}{2m},
$$

the allowed energy levels are

$$
\boxed{
E_n=\frac{n^2\pi^2\hbar^2}{2mL^2}
}
$$

where \(n\) is a positive integer.

This expression demonstrates that the energy is quantised and increases with the square of the quantum number.

For the \(n\)-th energy level,

$$
E_n=n^2E_1.
$$

Consequently,

$$
\frac{E_2}{E_1}=4
$$

and

$$
\frac{E_3}{E_1}=9.
$$

The energy levels also depend inversely on the square of the box length:

$$
E_n\propto\frac{1}{L^2}.
$$

Therefore, doubling the box length reduces every energy level to one-quarter of its original value.

### 2.4 Analytical Wavefunctions

The wavefunction must satisfy the normalisation condition

$$
\int_0^L|\psi_n(x)|^2\,dx=1.
$$

Applying this condition gives the normalised analytical wavefunctions

$$
\boxed{
\psi_n(x)=
\sqrt{\frac{2}{L}}
\sin\left(\frac{n\pi x}{L}\right)
}
$$

for \(0<x<L\).

The corresponding probability density is

$$
P_n(x)=|\psi_n(x)|^2.
$$

For the real-valued wavefunctions in this model,

$$
P_n(x)=\psi_n(x)^2.
$$

The probability density is non-negative and represents the spatial distribution of the probability of finding the particle.

The \(n\)-th wavefunction has \(n-1\) internal nodes. These are positions inside the box at which the wavefunction vanishes.

The ground state has no internal nodes, while the first excited state has one internal node at the centre of the box.

---

## 3. Numerical Methodology

### 3.1 Spatial Discretisation

A computer cannot directly represent a continuous wavefunction at infinitely many positions. Instead, the spatial domain is represented by a finite set of grid points.

For a box of length \(L\), the domain is divided into \(N\) equal intervals.

The grid spacing is

$$
\boxed{
h=\frac{L}{N}
}
$$

The grid points are defined by

$$
x_i=ih,
\qquad i=0,1,\ldots,N.
$$

Including both boundaries, the grid contains \(N+1\) points.

Since the boundary values are known to be zero, only the \(N-1\) interior wavefunction values need to be determined.

For the main calculation, \(N=200\). Therefore, the grid contains 201 points, including the boundaries, and 199 interior points.

### 3.2 Finite-Difference Approximation

The second derivative of the wavefunction can be approximated using the central finite-difference formula:

$$
\boxed{
\frac{d^2\psi}{dx^2}\bigg|_{x_i}
\approx
\frac{
\psi_{i+1}-2\psi_i+\psi_{i-1}
}{h^2}
}
$$

Substituting this approximation into the Schrödinger equation gives

$$
-\frac{\hbar^2}{2mh^2}
\left(
\psi_{i+1}-2\psi_i+\psi_{i-1}
\right)
=
E\psi_i.
$$

Define the constant

$$
C=\frac{\hbar^2}{2mh^2}.
$$

The equation can then be written as

$$
-C\psi_{i-1}+2C\psi_i-C\psi_{i+1}
=
E\psi_i.
$$

This equation connects the wavefunction at a grid point with its neighbouring values.

### 3.3 Hamiltonian Matrix Construction

The finite-difference equations can be combined into a matrix eigenvalue problem:

$$
\boxed{
\mathbf{H}\boldsymbol{\psi}_n
=
E_n\boldsymbol{\psi}_n
}
$$

For a system with three interior grid points, the Hamiltonian has the form

$$
\mathbf{H}
=
\begin{pmatrix}
2C & -C & 0\\
-C & 2C & -C\\
0 & -C & 2C
\end{pmatrix}.
$$

The diagonal elements are

$$
H_{ii}=\frac{\hbar^2}{mh^2},
$$

while the neighbouring off-diagonal elements are

$$
H_{i,i+1}=H_{i+1,i}
=
-\frac{\hbar^2}{2mh^2}.
$$

The resulting Hamiltonian is symmetric and tridiagonal.

Its eigenvalues approximate the allowed energy levels, and its eigenvectors represent the numerical wavefunctions at the interior grid points.

For \(N=200\), the Hamiltonian matrix has dimensions

$$
199\times199.
$$

### 3.4 Numerical Implementation

The calculations are implemented in Python using the following libraries:

- **NumPy:** Array operations and numerical calculations.
- **SciPy:** Eigenvalue calculations using a specialised tridiagonal matrix solver.
- **Matplotlib:** Visualisation of energy levels, wavefunctions, probability densities and convergence.

The analytical solution and numerical solution are calculated independently.

The numerical energy eigenvalues are subsequently compared with the analytical values to evaluate the accuracy of the finite-difference approximation.

---

## 4. Computational Parameters

The following physical constants and parameters are used in the calculations.

| Parameter | Symbol | Value |
|---|---|---|
| Reduced Planck constant | \(\hbar\) | \(1.054571817\times10^{-34}\,\mathrm{J\,s}\) |
| Electron mass | \(m\) | \(9.1093837139\times10^{-31}\,\mathrm{kg}\) |
| Box length | \(L\) | \(1.0\times10^{-9}\,\mathrm{m}\) |
| Electronvolt conversion | \(1\,\mathrm{eV}\) | \(1.602176634\times10^{-19}\,\mathrm{J}\) |
| Number of grid intervals | \(N\) | 200 |
| Number of interior grid points | \(N-1\) | 199 |

The numerical energy values are initially calculated in joules and then converted to electronvolts using

$$
E(\mathrm{eV})
=
\frac{E(\mathrm{J})}
{1.602176634\times10^{-19}}.
$$

---

## 5. Results and Discussion

### 5.1 Energy Eigenvalues

For an electron confined to a box of length \(1.0\,\mathrm{nm}\), the analytical ground-state energy is approximately

$$
E_1^{\mathrm{exact}}
\approx0.376030\,\mathrm{eV}.
$$

The numerical ground-state energy obtained using 200 grid intervals is approximately

$$
E_1^{\mathrm{num}}
\approx0.376022\,\mathrm{eV}.
$$

The small difference between these values arises from the finite-difference approximation.

The following table records the approximate numerical energy levels observed during development.

| Quantum number \(n\) | Numerical energy (eV) |
|---:|---:|
| 1 | 0.376022 |
| 2 | 1.503997 |
| 3 | 3.383645 |
| 4 | 6.014504 |
| 5 | 9.395923 |

These values should be checked against the output of the final analysis script before formal submission.

The numerical and analytical energy levels can be compared using the plot generated by the program:

`results/energy_comparison.png`

### 5.2 Energy Quantisation

The analytical energy levels satisfy

$$
E_n=n^2E_1.
$$

Therefore, the energy increases quadratically with the quantum number.

The increasing separation between successive energy levels is a characteristic feature of the infinite potential well.

The numerical calculation reproduces this general pattern.

### 5.3 Wavefunctions and Nodes

The analytical wavefunctions are sinusoidal functions satisfying the boundary conditions at both ends of the box.

For the \(n\)-th state, the number of internal nodes is

$$
N_{\mathrm{nodes}}=n-1.
$$

The ground-state wavefunction has no internal nodes. The first excited state has one node, while the second excited state has two.

The numerical wavefunctions are compared with their analytical counterparts using

`results/wavefunctions.png`

For this particular problem, the discrete eigenvectors of the uniform-grid finite-difference Hamiltonian have the same sampled sine-wave structure as the analytical solutions.

After appropriate normalisation and sign alignment, the numerical and analytical wavefunctions can agree at the grid points to floating-point precision.

This behaviour is a special mathematical property of the infinite potential well on a uniform grid. It should not be assumed for arbitrary potential-energy functions.

### 5.4 Probability Densities

The probability density is defined as

$$
P_n(x)=|\psi_n(x)|^2.
$$

It is non-negative and satisfies

$$
\int_0^L P_n(x)\,dx=1
$$

for a normalised state.

The probability density is zero wherever the wavefunction vanishes.

Consequently, the internal nodes of excited-state wavefunctions correspond to positions with zero probability density.

The calculated probability densities can be inspected using

`results/probability_densities.png`

### 5.5 Wavefunction Normalisation

For a continuous wavefunction, normalisation requires

$$
\int_0^L|\psi_n(x)|^2\,dx=1.
$$

The numerical eigensolver returns eigenvectors normalised according to the discrete condition

$$
\sum_{i=1}^{N-1}|v_i|^2=1.
$$

To obtain a continuum-scaled representation, the discrete eigenvectors are divided by the square root of the grid spacing:

$$
\psi(x_i)\approx\frac{v_i}{\sqrt{h}}.
$$

The normalisation integral can then be approximated numerically.

The resulting integral should be close to one, with any small discrepancy depending on the integration approximation and grid spacing.

---

## 6. Error Analysis and Convergence

### 6.1 Energy Error

The absolute energy error is defined as

$$
\Delta E_n
=
\left|
E_n^{\mathrm{num}}-E_n^{\mathrm{exact}}
\right|.
$$

The relative percentage error is

$$
\boxed{
\epsilon_n(\%)
=
\frac{
\left|E_n^{\mathrm{num}}-E_n^{\mathrm{exact}}\right|
}{
\left|E_n^{\mathrm{exact}}\right|
}
\times100
}
$$

This quantity provides a dimensionless measure of the difference between the numerical and analytical energy values.

### 6.2 Grid Convergence

The finite-difference method introduces a discretisation error because derivatives are approximated using values at discrete grid points.

For the central finite-difference approximation, the leading truncation error is proportional to \(h^2\).

At fixed box length,

$$
h=\frac{L}{N},
$$

so the leading energy error is expected to scale approximately as

$$
\boxed{
\epsilon\propto h^2\propto\frac{1}{N^2}
}
$$

when the calculation is in the asymptotic convergence regime.

Increasing the number of grid intervals reduces the spacing between grid points and generally improves the energy estimates.

However, a finer grid also increases the number of unknowns and the computational work.

The ground-state convergence results are saved by the program and visualised using

`results/energy_convergence.png`

The numerical data are stored in

`results/convergence_data.csv`

The convergence analysis provides an important test of the numerical method because it examines how the calculated energy changes as the spatial resolution is improved.

### 6.3 Energy Accuracy Versus Wavefunction Accuracy

Energy accuracy and wavefunction accuracy are related but distinct.

A small energy error does not, by itself, guarantee that a numerical wavefunction accurately represents the physical state.

For the infinite potential well on a uniform grid, the numerical eigenvectors reproduce the analytical sine-wave values at the sampled positions exceptionally closely.

Therefore, the sampled wavefunction comparison is unusually favourable for this system.

For more complicated potentials, numerical wavefunctions generally require separate and careful validation.

---

## 7. Limitations

Although the infinite potential well is useful for studying computational quantum mechanics, the model has several limitations.

1. **Idealised confinement:** The potential barriers are assumed to be infinitely high, which is an idealisation.
2. **One-dimensional system:** The calculation does not account for motion in three spatial dimensions.
3. **Finite-difference error:** The numerical Hamiltonian approximates the continuous differential operator, introducing discretisation error.
4. **Grid dependence:** Numerical accuracy depends on the selected grid resolution.
5. **Special wavefunction agreement:** The close agreement between numerical and analytical wavefunctions at grid points is specific to this model and should not be generalised.
6. **Limited chemical realism:** The infinite potential well is a foundational model rather than a realistic molecular electronic-structure calculation.

Despite these limitations, the model provides a useful introduction to numerical eigenvalue problems and computational quantum mechanics.

---

## 8. Conclusion

This project investigated the one-dimensional infinite potential well using analytical quantum mechanics and numerical computation.

The time-independent Schrödinger equation was discretised using a central finite-difference approximation, producing a symmetric tridiagonal Hamiltonian matrix.

Numerical diagonalisation yielded energy eigenvalues and eigenvectors that could be compared with the analytical solutions.

The calculated energy levels reproduced the expected quantisation pattern and agreed closely with the analytical results. Wavefunction and probability-density analysis illustrated the relationship between quantum states, internal nodes and spatial probability distributions.

The grid-convergence analysis demonstrated the importance of spatial resolution in numerical calculations.

Overall, the project established a computational workflow involving differential-equation discretisation, matrix construction, eigenvalue calculations, visualisation and numerical validation.

These techniques provide a foundation for studying more complex quantum systems, including the quantum harmonic oscillator, hydrogen atom and molecular electronic-structure models.

---

## 9. Reproducibility

All calculations are implemented in Python.

The required packages are listed in `requirements.txt`.

Activate the virtual environment and install the dependencies using

```bash
python -m pip install -r requirements.txt
```

Run the scripts from the project root directory:

```bash
python src/analytical_solution.py
python src/numerical_solution.py
python src/wavefunction_validation.py
python src/final_analysis.py
```

The scripts generate numerical results, plots and CSV files in the `results/` directory.

For reproducibility, the final project should preserve the code, computational parameters, generated results and documentation together.

The Python, NumPy, SciPy and Matplotlib versions used for the final calculations should also be recorded.

---

## 10. Future Work

The methods developed in this project can be extended to more complicated quantum-mechanical systems.

Potential future investigations include:

- The quantum harmonic oscillator and its vibrational energy levels.
- Finite potential wells and quantum tunnelling.
- The hydrogen atom and hydrogenic orbitals.
- More general potential-energy functions.
- Improved numerical methods for solving the Schrödinger equation.
- Applications to molecular vibrations and electronic structure.

These extensions will allow the computational framework developed here to be applied to increasingly realistic physical and chemical systems.