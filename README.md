<div align="center">

# 🔬 PHY341/745 — Physics Computer Lab

[![NISER](https://img.shields.io/badge/NISER-Physics-blue?style=for-the-badge)](https://www.niser.ac.in)
[![Python](https://img.shields.io/badge/Python-3.x-3776AB?style=for-the-badge&logo=python&logoColor=white)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow?style=for-the-badge)](LICENSE)

**Numerical Methods in Physics — implemented from scratch**

*Instructor: Prof. Subhasish Basak · Student: Omm Biswal (Roll No. 2411122)*

</div>

---

## 📋 Table of Contents

- [Course Overview](#-course-overview)
- [Marks & Rules](#-marks--rules)
- [Programming Rules](#-programming-rules)
- [Repository Structure](#-repository-structure)
- [Topics](#-topics)
  - [1 · Random Number Generation](#1--random-number-generation)
  - [2 · Gauss-Jordan Elimination](#2--gauss-jordan-elimination)
  - [3 · LU Decomposition](#3--lu-decomposition)
  - [4 · Iterative Methods — Jacobi & Gauss-Seidel](#4--iterative-methods--jacobi--gauss-seidel)
  - [5 · Root Finding](#5--root-finding)
  - [6 · Numerical Integration](#6--numerical-integration)
  - [7 · Differential Equations (ODEs & PDEs)](#7--differential-equations-odes--pdes)
  - [8 · Eigenvalues & Eigenvectors](#8--eigenvalues--eigenvectors)
  - [9 · Interpolation & Data Fitting](#9--interpolation--data-fitting)
- [Reference Book](#-reference-book)

---

## 🎓 Course Overview

This lab course covers the essential **numerical methods** used in computational physics — all implemented from scratch in Python or C/C++, no solver libraries allowed for the main algorithms.

| | |
|---|---|
| **Credits** | 4 (Int-MSc Sem-V) / 4 (Int-MSc-PhD Sem-III) |
| **Schedule** | LH-1, Wednesday & Thursday @ 2:30 PM |
| **Instructor** | Prof. Subhasish Basak · SPS-507 · sbasak@niser.ac.in |
| **TAs** | Hemant Lohumi (SPS-503B) · Chinmoy Samanta (SPS-412B) |
| **Platform** | Jupyter / VS Code · GitHub recommended |
| **Languages** | Python or C/C++ only |

---

## 📊 Marks & Rules

| Component | Marks |
|-----------|-------|
| Assignments | 35 |
| Midsem | 15 |
| DIY Project | 20 |
| Endsem | 30 |

- Assignments close at deadline — no late submissions accepted.
- Suspected copying → zero for **all** parties, no chance to defend.
- Midsem is the first class after Midsem recess (not a formal lab exam day).
- Endsem is full-syllabus and **must use your own library functions** — no system built-ins for major algorithms.
- No more than 10 % of the class receives AA.

---

## 💻 Programming Rules

1. Every file has a **header comment**: purpose, name, roll number.
2. **No hardcoded or interactive input** — read everything from external files (matrix size, precision, filenames, etc.).
3. Output goes to a **separate file** that must exactly match what the code produces. Comments may be appended below, clearly separated.
4. Code must have **reasonable comments** — enough for someone else to follow, not a wall of text.
5. **Loop-based initialisation** — no inline array setup for things like reading matrices.
6. **One growing library file** (`MyLib.py` / `mylib.c`) holds all functions: PRNG, Gauss-Jordan, LU, root finders, integrators, etc. Never scatter them into main files.
7. The **main file is thin**: only I/O, memory, file handles, and calls into the library.
8. **No system built-ins** for major algorithms — `random`, `numpy.linalg.solve`, `scipy.integrate` are banned. Basic math (`sqrt`, `sin`, `log`) is fine.

---

## 📁 Repository Structure

```
CS-341-Computational-Lab/
│
├── CS-341-slides/          # Lecture slides (PDFs)
│   ├── intro.pdf           # Course intro & programming practices
│   ├── random.pdf          # Random number generation
│   ├── gaussj.pdf          # Gauss-Jordan elimination
│   ├── ludecomp.pdf        # LU decomposition & Cholesky
│   ├── itermat.pdf         # Jacobi, Gauss-Seidel, SOR
│   ├── root.pdf            # Root finding methods
│   ├── integ.pdf           # Numerical integration
│   ├── diffeq.pdf          # ODEs, BVPs, and PDEs (Heat & ADR equations)
│   ├── eigv.pdf            # Eigenvalues & Eigenvectors (Power Iteration & Deflation)
│   └── datafit.pdf         # Interpolation & Data Fitting (Lagrange & Least Squares)
│
├── MyLib.py                # Core library: PRNG, MyComplex
├── MyLib_Advanced.py       # Extended library: Solve, Newton-Raphson,
│                           #   Fixed-Point, Deriv, GaussPRNG
├── PRNG.py                 # Standalone PRNG module
├── PRNG_Matrices.py        # PRNG + matrix utilities
│
├── Daily_Work/
│   ├── Day_0/ … Day_19/    # Daily class work, scripts, and outputs
│   └── (Day_19/MyLib.py)   # Most complete library snapshot
│
└── README.md
```

---

## 📐 Topics

> **Note on equations:** GitHub renders LaTeX in Markdown via `$...$` (inline) and `$$...$$` (block). All formulas below use this syntax.

---

### 1 · Random Number Generation

> 📄 [`CS-341-slides/random.pdf`](CS-341-slides/random.pdf)

**Goal:** Write your own pseudo-random number generator and use it for every assignment, exam, and DIY — never call `random` or `numpy.random`.

#### True vs Pseudo RNG

| | True RNG | Pseudo RNG |
|---|---|---|
| Source | Physical phenomena (radioactivity, thermal noise…) | Deterministic algorithm |
| Speed | Slow | Extremely fast |
| Portability | Low | High |
| Period | Infinite | Finite but tunably long |

#### Logistic Map (quick & dirty)

$$x_{i+1} = c \, x_i \,(1 - x_i)$$

Bad for most $c$; shows no discernible pattern only near $c = 3.98,\ x_0 = 0.1$.

#### Linear Congruential Generator (LCG)

$$x_{i+1} = (a\,x_i + c)\bmod m, \qquad r_i = x_i / m \in [0,1)$$

| Parameter | Numerical Recipes | gcc |
|---|---|---|
| $a$ | 1 664 525 | 1 103 515 245 |
| $c$ | 1 013 904 223 | 12 345 |
| $m$ | $2^{32}$ | $2^{31}$ |

**Hull-Dobell Theorem** — LCG has full period $m$ iff $c \neq 0$ and:
1. $c$ is coprime to $m$
2. $a - 1$ is divisible by every prime factor of $m$
3. $a - 1$ is divisible by 4 if $m$ is divisible by 4

#### Non-uniform distributions

| Distribution | Formula |
|---|---|
| Uniform $[a,b)$ | $u = a + (b-a)\,x$ |
| Exponential with rate $\alpha$ | $y = -\tfrac{1}{\alpha}\ln x$ |
| Gaussian (Box-Muller) | $z = \sqrt{-2\ln u_1}\cos(2\pi u_2)$ |

#### Testing randomness

- **Basic:** correlations and moments
- **Advanced:** chi-square, Kolmogorov-Smirnov
- Connected correlation: $\varepsilon(n,N) = \frac{1}{N}\sum x_i x_{i+n} - \left(\frac{1}{N}\sum x_i\right)\!\left(\frac{1}{N}\sum x_{i+n}\right) \to 0$

#### Applications

- **Estimating $\pi$:** throw random $(x,y) \in [0,1]^2$; fraction with $x^2+y^2 \le 1$ approximates $\pi/4$.
- **Radioactive decay:** draw $r \in [0,1]$; if $r < P_A = \lambda_A \,\Delta t$, decay $A \to B$.

#### Code — `MyLib.py`

```python
def PRNG(seed=42):          # Linear Congruential Generator
    a = 1664525.0
    c = 1013904223.0
    m = 2**32
    seed = (a * seed + c) % m
    return float(seed / m)
```

```python
# MyLib_Advanced.py — Gaussian RNG via Box-Muller
def GaussPRNG(seed=42.0, iter=1):
    u = list(PRNG(seed, 2 * iter))
    for i in range(iter):
        u1, u2 = u[2*i], u[2*i+1]
        yield math.sqrt(-2 * math.log(u1)) * math.cos(2 * math.pi * u2)
```

---

### 2 · Gauss-Jordan Elimination

> 📄 [`CS-341-slides/gaussj.pdf`](CS-341-slides/gaussj.pdf)

**Use case:** Small, dense matrices. Memory: $O(N^2)$. Time: $O(N^3)$.

#### Augmented matrix → RREF

Write $[A \mid \mathbf{b}]$ and reduce to Reduced Row Echelon Form $[I \mid \tilde{\mathbf{b}}]$ — reading off $x_i = \tilde{b}_i$.

RREF conditions:
1. Zero rows at the bottom.
2. Each pivot is strictly to the right of the one above it.
3. Each pivot is exactly 1.
4. All other entries in a pivot's column are 0.

#### Partial pivoting (essential)

Before reducing column $j$, swap the current row with whichever row below has the **largest absolute value** in column $j$. This prevents division by near-zero values and keeps rounding errors small. Track the number of swaps $n$ for the determinant.

$$\det(A) = (-1)^n \prod a'_{ii}$$

#### Matrix inversion

Augment with the identity matrix: $[A \mid I] \xrightarrow{\text{G-J}} [I \mid A^{-1}]$.

> ⚠ For $N \gtrsim 10$ rounding errors build up fast. Use LU decomposition instead.

---

### 3 · LU Decomposition

> 📄 [`CS-341-slides/ludecomp.pdf`](CS-341-slides/ludecomp.pdf)

Same $O(N^3)$ cost as Gauss-Jordan, but reusable for multiple right-hand sides and gives the determinant for free.

#### Doolittle factorisation ($\ell_{ii} = 1$)

Factorise $A = LU$ in-place (overwrite $a_{ij}$ with $u_{ij}$ for $i \le j$, $\ell_{ij}$ for $i > j$):

$$u_{ij} = a_{ij} - \sum_{k=0}^{i-1} \ell_{ik}\,u_{kj}, \quad i \le j$$

$$\ell_{ij} = \frac{1}{u_{jj}}\!\left(a_{ij} - \sum_{k=0}^{j-1} \ell_{ik}\,u_{kj}\right), \quad i > j$$

#### Solving $A\mathbf{x} = \mathbf{b}$

Split into two triangular solves:

**1. Forward substitution** — solve $L\mathbf{y} = \mathbf{b}$:
$$y_i = b_i - \sum_{j=0}^{i-1} \ell_{ij}\,y_j$$

**2. Backward substitution** — solve $U\mathbf{x} = \mathbf{y}$:
$$x_i = \frac{1}{u_{ii}}\!\left(y_i - \sum_{j=i+1}^{N-1} u_{ij}\,x_j\right)$$

#### Determinant & inverse for free

$$\det(A) = (-1)^n \prod_{i} u_{ii}$$

For the inverse, run forward + backward substitution over each column of $I$.

#### Cholesky decomposition

For **symmetric positive-definite** matrices: $A = LL^T$

$$\ell_{ii} = \sqrt{a_{ii} - \sum_{j < i} \ell_{ij}^2}, \qquad \ell_{ij} = \frac{1}{\ell_{ii}}\!\left(a_{ij} - \sum_{k=0}^{i-1} \ell_{ik}\,\ell_{jk}\right),\ j > i$$

About **twice as fast** as plain LU. Used in covariance matrix decomposition and Monte Carlo simulation.

#### Code — `MyLib_Advanced.py`

```python
def Solve(A, b):
    """Auto-detects symmetric matrices and uses Cholesky; falls back to LU."""
    A = Copy_Mat(A)
    n = len(A)
    symmetric = all(A[i][j] == A[j][i] for i in range(n) for j in range(n))
    if symmetric:
        L  = Cholesky(A)
        y  = forward(L, b)
        Lt = Cholesky_read(Cholesky_compress(L), 1)
        x  = backward(Lt, y)
    else:
        LUmat, perm = LU(A)
        bperm = [b[p] for p in perm]
        y = forward(LUmat, bperm)
        x = backward(LUmat, y)
    return x
```

---

### 4 · Iterative Methods — Jacobi & Gauss-Seidel

> 📄 [`CS-341-slides/itermat.pdf`](CS-341-slides/itermat.pdf)

**Use case:** Large or sparse matrices where $O(N^3)$ direct methods are too slow or memory-hungry.

**Convergence is guaranteed** when $A$ is strictly diagonally dominant:
$$|a_{ii}| > \sum_{j \neq i} |a_{ij}|$$

#### Jacobi method

Decompose $A = D + (L+U)$. Iterate using **only old values**:

$$x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j \neq i} a_{ij}\,x_j^{(k)} \right)$$

- Needs **two storage vectors** (cannot overwrite $x^{(k)}$ with $x^{(k+1)}$ mid-sweep).
- Main advantage: trivially **parallelisable** — great for large sparse systems.

#### Gauss-Seidel method

Uses newly computed values immediately:

$$x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j < i} a_{ij}\,x_j^{(k+1)} - \sum_{j > i} a_{ij}\,x_j^{(k)} \right)$$

- Only **one storage vector** needed.
- Converges roughly **2× faster** than Jacobi in practice.

#### Successive Over-Relaxation (SOR)

Blend Gauss-Seidel with the previous iterate using weight $\omega$:

$$x_i^{(k+1)} = (1-\omega)\,x_i^{(k)} + \frac{\omega}{a_{ii}} \left( b_i - \sum_{j<i} a_{ij}\,x_j^{(k+1)} - \sum_{j>i} a_{ij}\,x_j^{(k)} \right)$$

| $\omega$ | Effect |
|---|---|
| $= 1$ | Normal Gauss-Seidel |
| $1 < \omega < 2$ | Over-relaxation — speeds up slow convergence |
| $0 < \omega < 1$ | Under-relaxation — stabilises oscillating solutions |

> **Example:** For a tridiagonal $n \times n$ system ($n \ge 3$), Gauss-Seidel takes ~150 iterations; SOR with $\omega \approx 1.57$ takes ~25.

**Convergence check** for all iterative methods:
$$\|\mathbf{x}^{(k+1)} - \mathbf{x}^{(k)}\| < \varepsilon$$

---

### 5 · Root Finding

> 📄 [`CS-341-slides/root.pdf`](CS-341-slides/root.pdf)

**Goal:** Solve $f(x) = 0$ numerically.

**Bracketing:** Find $[a, b]$ where $f(a) \cdot f(b) < 0$. By the Intermediate Value Theorem, at least one root lies inside.

#### Bisection method

Simplest; guaranteed to converge if properly bracketed. Halves the interval each step.

$$c = \frac{a+b}{2}, \quad \text{then replace } a \text{ or } b \text{ with } c \text{ depending on sign of } f(c)$$

After $n$ steps: $|b_n - a_n| = |b-a|/2^n \le \varepsilon$. Slow — linear convergence.

#### Regula Falsi (False Position)

Draws a secant through $(a, f(a))$ and $(b, f(b))$; the x-intercept $c$ is the new bracket point:

$$c_n = b_n - \frac{(b_n - a_n)\,f(b_n)}{f(b_n) - f(a_n)}$$

Always converges; faster than bisection, especially for nearly-linear functions.

#### Newton-Raphson method

Requires $f'(x)$; no bracketing needed. **Quadratic convergence** near the root (significant digits roughly double each step):

$$x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$$

**Multivariable Newton-Raphson** — uses the Jacobian $J_{pq} = \partial f_p / \partial x_q$:

$$\mathbf{x}^{(i+1)} = \mathbf{x}^{(i)} - J^{-1}(\mathbf{x}^{(i)})\,\mathbf{f}(\mathbf{x}^{(i)})$$

#### Fixed-Point / Picard iteration

Rewrite $f(x) = 0$ as $x = g(x)$. Iterate $x_{n+1} = g(x_n)$.

Converges (stable fixed point) when $|g'(x_0)| < 1$. Multiple valid choices of $g$ exist for one $f$ — some converge, some don't.

#### Laguerre's method (polynomials)

For a degree-$n$ polynomial $P(x)$. Two-stage process:

1. **Laguerre iteration** — find root $\alpha_1$:

$$G = \frac{P'(\beta_k)}{P(\beta_k)}, \quad H = G^2 - \frac{P''(\beta_k)}{P(\beta_k)}, \quad a = \frac{n}{G \pm \sqrt{(n-1)(nH - G^2)}}$$

Choose the sign giving the **larger** denominator. Set $\beta_{k+1} = \beta_k - a$.

2. **Deflation** (synthetic division): divide out $(x - \alpha_1)$ to get a degree $n-1$ polynomial, then repeat.

#### Code — `MyLib_Advanced.py`

```python
def Deriv(f, x, i=None, h=1e-6):
    """Finite-difference derivative estimate (secant method)."""
    if type(x) is not list:          # 1D
        return (f(x + h) - f(x)) / h
    else:                            # multivariable: partial w.r.t. x[i]
        x_fwd = x[:]
        x_fwd[i] += h
        return (f(x_fwd) - f(x)) / h

def Newton_Raphson(F, x0, h=1e-6, tolerance=1e-6, max_iter=10000):
    """
    Derivative-free Newton-Raphson (uses Deriv() internally).
    F  : scalar function, or list of functions for multivariable case.
    x0 : float (1D) or list of floats (multivariable).
    """
    iters = 0
    if type(x0) is not list:        # ── 1D ──
        x = x0
        while iters < max_iter:
            x_new = x - F(x) / Deriv(F, x, h=h)
            iters += 1
            if abs(x_new - x) < tolerance:
                break
            x = x_new
    else:                           # ── multivariable ──
        n, x = len(x0), x0[:]
        while iters < max_iter:
            fx  = [F[i](x) for i in range(n)]
            Jx  = [[Deriv(F[i], x, j, h) for j in range(n)] for i in range(n)]
            Jinv = invert(Jx)
            dx  = [sum(Jinv[i][j] * fx[j] for j in range(n)) for i in range(n)]
            x_new = [x[i] - dx[i] for i in range(n)]
            iters += 1
            if math.sqrt(sum((x_new[i]-x[i])**2 for i in range(n))) < tolerance:
                break
            x = x_new
    return x_new
```

---

### 6 · Numerical Integration

> 📄 [`CS-341-slides/integ.pdf`](CS-341-slides/integ.pdf)

Approximate $\displaystyle\mathcal{I} = \int_a^b f(x)\,dx \approx \sum_{n=1}^{N} w(x_n)\,f(x_n)$

with $N$ evaluation points and weights $w(x_n)$ derived from Lagrange interpolation.

#### Midpoint rule

Evaluate $f$ at the **midpoint** of each sub-interval of width $h = (b-a)/N$:

$$\mathcal{M}_N = h \sum_{n=1}^{N} f\!\left(x_n\right), \quad x_n = a + \left(n - \tfrac{1}{2}\right)h$$

Error bound: $\displaystyle\lvert\mathcal{M}_N - \mathcal{I}\rvert \le \frac{(b-a)^3}{24N^2}\,\lvert f''\rvert_{\max}$

#### Trapezoidal rule

Connect consecutive $f(x_n)$ with straight lines:

$$\mathcal{T}_N = \frac{h}{2}\!\left[f(x_0) + 2f(x_1) + \cdots + 2f(x_{N-1}) + f(x_N)\right]$$

Error bound: $\displaystyle\lvert\mathcal{T}_N - \mathcal{I}\rvert \le \frac{(b-a)^3}{12N^2}\,\lvert f''\rvert_{\max}$

> ⚠ Surprisingly, the trapezoidal rule can be *less* accurate than midpoint for strictly concave or convex functions.

#### Simpson's 1/3 rule

Uses parabolic interpolation over pairs of sub-intervals (**even** $N$ required):

$$\mathcal{S}_N = \frac{h}{3}\!\left[f(x_0) + 4f(x_1) + 2f(x_2) + 4f(x_3) + \cdots + 4f(x_{N-1}) + f(x_N)\right]$$

Weights follow the pattern $1, 4, 2, 4, 2, \ldots, 4, 1$.

Error bound: $\displaystyle\lvert\mathcal{S}_N - \mathcal{I}\rvert \le \frac{(b-a)^5}{180N^4}\,\lvert f^{(4)}\rvert_{\max}$

Useful identity: $\mathcal{S}_{2N} = \tfrac{2}{3}\mathcal{M}_N + \tfrac{1}{3}\mathcal{T}_N$

#### Error summary

| Method | Error scales as | Requires |
|---|---|---|
| Midpoint | $O(N^{-2})$ | $f''$ bounded |
| Trapezoidal | $O(N^{-2})$ | $f''$ bounded |
| Simpson | $O(N^{-4})$ | $f''''$ bounded, even $N$ |
| Monte Carlo | $O(N^{-1/2})$ | uniform pRNG |

#### Monte Carlo integration

For uniform draws $X_i \in [a,b]$ via your pRNG ($X = a + (b-a)\xi,\ \xi\in[0,1]$):

$$\mathcal{F}_N = \frac{b-a}{N}\sum_{i=1}^N f(X_i), \qquad \sigma_f^2 = \frac{1}{N}\sum f(X_i)^2 - \left(\frac{1}{N}\sum f(X_i)\right)^{\!2}$$

Converges slowly ($\sim 1/\sqrt{N}$) but is the **only practical method** for high-dimensional integrals.

#### Gaussian quadrature

Uses $n$ non-equally-spaced nodes (roots of orthogonal polynomials) + optimal weights — exact for all polynomials up to degree $2n-1$.

| Integral type | Polynomial |
|---|---|
| $\int_{-1}^{1}$ | Legendre $P_n(x)$ |
| $\int_0^\infty$ | Laguerre $L_n(x)$ |
| $\int_{-\infty}^{\infty}$ | Hermite $H_n(x)$ |

For $\int_a^b$, change variables: $t = \tfrac{b-a}{2}x + \tfrac{b+a}{2}$.

**Efficiency comparison** for $\int_{-1}^{1} x\,e^x\,dx = 2/e \approx 0.735\,758\,88$:

| Simpson $N$ | Result | Gauss $n$ | Result |
|---|---|---|---|
| 40 | 0.735 759 23 | 4 | 0.735 756 50 |
| 80 | 0.735 758 90 | 5 | 0.735 758 87 |
| 120 | 0.735 758 88 | **6** | **0.735 758 88** ✓ |

Gaussian quadrature matches Simpson at $N=120$ using only **6 points**.

---

### 7 · Differential Equations (ODEs & PDEs)

> 📄 [`CS-341-slides/diffeq.pdf`](CS-341-slides/diffeq.pdf)

Covers initial value problems (IVPs) for ordinary differential equations, boundary value problems (BVPs) via shooting and finite difference schemes, and partial differential equations (PDEs) focusing on 1D diffusion/heat equations and Crank-Nicolson schemes.

#### Initial Value Problems (IVPs)

First-order ODE: $\displaystyle \frac{dy}{dx} = f(x, y), \quad y(x_0) = y_0$

Higher-order ODEs are converted into coupled systems of first-order ODEs. For example, a damped harmonic oscillator $\ddot{x} + \mu \dot{x} + \omega^2 x = 0$:
$$\frac{dx}{dt} = v, \qquad \frac{dv}{dt} = -\mu v - \omega^2 x, \quad x(t_0) = x_0, \ v(t_0) = v_0$$

#### Forward (Explicit) Euler

Evaluates the derivative at the beginning of the step:
$$y_{n+1} = y_n + h f(x_n, y_n) + O(h^2)$$

- **Local truncation error:** $O(h^2)$; **global accumulation error:** $O(h)$.
- **Stability issue:** Strictly conditionally stable. The solution can rapidly diverge (blow up) unless step size $h$ is sufficiently small.

#### Backward (Implicit) Euler

Evaluates the derivative at the *end* of the step:
$$y_{n+1} = y_n + h f(x_{n+1}, y_{n+1})$$

- Implicit algebraic equation in $y_{n+1}$ — solved via root-finding (Newton-Raphson) at every step:
  $$G(z) = z - y_n - h f(x_{n+1}, z) = 0$$
- **Unconditionally stable:** Allows significantly larger step sizes without numerical divergence, compensating for the extra computational cost per step.

#### Predictor-Corrector Method (Heun's)

Combines forward Euler as a *predictor* with the trapezoidal rule as a *corrector*:
1. **Predict:** $y^p_{n+1} = y_n + h f(x_n, y_n) = y_n + k_1$
2. **Correct:** $y^c_{n+1} = y_n + \frac{h}{2}\!\left[f(x_n, y_n) + f(x_{n+1}, y^p_{n+1})\right] = y_n + \frac{1}{2}(k_1 + k_2)$
- Local error: $O(h^3)$. Improves accuracy and stability over plain forward Euler.

#### Runge-Kutta 4th Order (RK4)

The primary workhorse for ODE initial value problems. Approximates the integral $\int_{x_n}^{x_n+h} f(x, y)\,dx$ using a Simpson-like weighted average of four trial slopes:

$$k_1 = h f(x_n, y_n)$$
$$k_2 = h f\!\left(x_n + \frac{h}{2},\, y_n + \frac{k_1}{2}\right)$$
$$k_3 = h f\!\left(x_n + \frac{h}{2},\, y_n + \frac{k_2}{2}\right)$$
$$k_4 = h f(x_n + h,\, y_n + k_3)$$

$$y_{n+1} = y_n + \frac{1}{6}\!\left(k_1 + 2k_2 + 2k_3 + k_4\right) + O(h^5)$$

- Global error is $O(h^4)$, enabling high accuracy even with relatively coarse steps.
- **Coupled Systems:** Extended directly to $m$-dimensional state vectors $\mathbf{Y} = [y_1, y_2, \ldots, y_m]^T$:
  $$\mathbf{Y}_{n+1} = \mathbf{Y}_n + \frac{1}{6}\!\left(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4\right)$$
  Used to simulate non-linear dynamical and chaotic systems, including the **Lorenz attractor** ($\dot{x} = \sigma(y-x), \dot{y} = x(\rho-z)-y, \dot{z} = xy-\beta z$) and **Rössler system**.

#### Boundary Value Problems (BVPs)

Second-order ODE with conditions specified at boundaries $x = a$ and $x = b$:
$$y'' = f(x, y, y'), \quad a \le x \le b$$
- **Dirichlet conditions:** $y(a) = \alpha, \quad y(b) = \beta$
- **Neumann conditions:** $y'(a) = Y_0, \quad y'(b) = Y_N$

##### 1. Shooting Method
Transforms the BVP into an initial value problem by guessing the initial slope $z(a) = y'(a) = \zeta$:
1. Solve the coupled system $\frac{dy}{dx} = z, \ \frac{dz}{dx} = f(x, y, z)$ using RK4 from $x = a$ to $x = b$.
2. Bracket the target boundary value $\beta$: find $\zeta_l$ (undershooting $y_{\zeta_l}(b) < \beta$) and $\zeta_h$ (overshooting $y_{\zeta_h}(b) > \beta$).
3. Interpolate the next slope via Lagrange interpolation:
   $$\zeta = \zeta_l + \frac{\zeta_h - \zeta_l}{y_{\zeta_h}(b) - y_{\zeta_l}(b)}\!\left(\beta - y_{\zeta_l}(b)\right)$$
4. Repeat until $|y_\zeta(b) - \beta| < \varepsilon$.

##### 2. Finite Difference / Finite Element Method (1D)
Discretizes the domain $x_n = a + n h$ with $h = (b-a)/N$. Approximates $y''(x_n)$ via central differences:
$$y''(x_n) \approx \frac{y_{n+1} - 2y_n + y_{n-1}}{h^2}$$
For interior nodes $n = 1, \dots, N-1$, this yields a tridiagonal, symmetric positive-definite linear system:
$$-y_{n-1} + 2y_n - y_{n+1} = -h^2 f(x_n)$$
with Dirichlet boundary conditions $\alpha, \beta$ shifted to the right-hand side. Solved directly using **Cholesky decomposition** ($O(N)$ operations).

#### Partial Differential Equations (PDEs)

General linear 2nd-order PDE in two variables $A u_{xx} + B u_{xy} + C u_{yy} = F$:
- **Elliptic ($B^2 - 4AC < 0$):** Steady-state equilibrium (Laplace $\nabla^2 u = 0$, Poisson $\nabla^2 u = \rho$).
- **Parabolic ($B^2 - 4AC = 0$):** Diffusion and dissipation ($u_t = D u_{xx}$).
- **Hyperbolic ($B^2 - 4AC > 0$):** Wave propagation ($u_{tt} = c^2 u_{xx}$).

##### 1D Heat / Diffusion Equation ($u_t = u_{xx}$)
Grid discretization: $x_i = i h_x$, $t_j = j h_t$, with parameter $\alpha = h_t / h_x^2$.

- **Explicit FTCS (Forward-Time Central-Space):**
  $$\frac{u_{i, j+1} - u_{i, j}}{h_t} = \frac{u_{i+1, j} - 2u_{i, j} + u_{i-1, j}}{h_x^2} \implies u_{i, j+1} = \alpha u_{i-1, j} + (1 - 2\alpha)u_{i, j} + \alpha u_{i+1, j}$$
  In vector form: $\mathbf{V}_{j+1} = A \mathbf{V}_j$.
  > ⚠ **Stability criterion:** Requires $\alpha = \frac{h_t}{h_x^2} \le \frac{1}{2}$. If $h_t > \frac{1}{2}h_x^2$, numerical oscillations grow exponentially.

- **Implicit BTCS (Backward-Time Central-Space):**
  Uses backward time derivative: $-\alpha u_{i-1, j} + (1 + 2\alpha)u_{i, j} - \alpha u_{i+1, j} = u_{i, j-1}$.
  In vector form: $A \mathbf{V}_j = \mathbf{V}_{j-1}$. Unconditionally stable for all $\alpha$.

- **Crank-Nicolson Scheme:**
  Averages central spatial differences between time steps $j$ and $j-1$:
  $$\frac{u_{i, j} - u_{i, j-1}}{h_t} = \frac{1}{2}\!\left[\frac{u_{i+1, j} - 2u_{i, j} + u_{i-1, j}}{h_x^2} + \frac{u_{i+1, j-1} - 2u_{i, j-1} + u_{i-1, j-1}}{h_x^2}\right]$$
  In matrix notation:
  $$(2I + \alpha B)\mathbf{V}_j = (2I - \alpha B)\mathbf{V}_{j-1}$$
  where $B$ is the tridiagonal matrix with main diagonal $2$ and off-diagonals $-1$.
  - **Properties:** Unconditionally stable, second-order accurate in both time and space ($O(h_t^2 + h_x^2)$).
  - Generalizes directly to **Advection-Diffusion-Reaction (ADR)** equations: $u_t = a u_{xx} + b u_x + c u$.

#### Code — `MyLib.py`

```python
def RK4(f, x0, Y0, x_end, h):
    """Runge-Kutta 4th order for 1st order ODE or coupled system dY/dx = f(x, Y)."""
    xs, ys = [x0], [Y0]
    x, Y = x0, Y0
    system = isinstance(Y0, list)
    while x < x_end - 1e-12:
        step = min(h, x_end - x)
        if system:
            n = len(Y)
            k1 = f(x, Y)
            k2 = f(x + step/2, [Y[i] + step/2 * k1[i] for i in range(n)])
            k3 = f(x + step/2, [Y[i] + step/2 * k2[i] for i in range(n)])
            k4 = f(x + step,   [Y[i] + step   * k3[i] for i in range(n)])
            Y = [Y[i] + (step/6) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(n)]
        else:
            k1 = f(x, Y)
            k2 = f(x + step/2, Y + step/2 * k1)
            k3 = f(x + step/2, Y + step/2 * k2)
            k4 = f(x + step,   Y + step   * k3)
            Y += (step/6) * (k1 + 2*k2 + 2*k3 + k4)
        x += step
        xs.append(x)
        ys.append(Y)
    return xs, ys
```

---

### 8 · Eigenvalues & Eigenvectors

> 📄 [`CS-341-slides/eigv.pdf`](CS-341-slides/eigv.pdf)

Addresses the algebraic eigenvalue problem $A\mathbf{x} = \lambda\mathbf{x}$.

#### Numerical Considerations

Analytically, eigenvalues are roots of the characteristic polynomial $p_A(\lambda) = \det(A - \lambda I) = 0$. In practice, polynomial root-finding for eigenvalues is **numerically unstable and ill-conditioned** (Wilkinson's phenomenon: tiny perturbations in matrix entries cause massive changes in polynomial roots). Iterative numerical algorithms exploiting matrix structure (symmetry, sparsity, tridiagonality) are essential.

#### Power Iteration Method

An iterative method designed to extract the **dominant eigenvalue** (largest absolute magnitude $|\lambda_1|$) and its corresponding **dominant eigenvector** $\mathbf{v}_1$.

**Requirements:**
1. Matrix $A$ has $n$ linearly independent eigenvectors $\{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$.
2. Eigenvalues can be ordered in magnitude: $|\lambda_1| > |\lambda_2| \ge \cdots \ge |\lambda_n|$.

**Algorithm:**
1. Choose an initial unit vector $\mathbf{x}_0$ with non-zero projection along $\mathbf{v}_1$ ($\mathbf{x}_0 = \sum c_i \mathbf{v}_i, \ c_1 \ne 0$).
2. Successive matrix-vector multiplications magnify the dominant component:
   $$A^k \mathbf{x}_0 = c_1 \lambda_1^k \mathbf{v}_1 + c_2 \lambda_2^k \mathbf{v}_2 + \cdots = c_1 \lambda_1^k \left[\mathbf{v}_1 + \sum_{i=2}^n \frac{c_i}{c_1}\left(\frac{\lambda_i}{\lambda_1}\right)^{\!k}\mathbf{v}_i\right] \approx c_1 \lambda_1^k \mathbf{v}_1$$
3. Estimate the dominant eigenvalue at each step using the **Rayleigh Quotient**:
   $$\lambda^{(k)} = \frac{A\mathbf{x}_k \cdot \mathbf{x}_k}{\mathbf{x}_k \cdot \mathbf{x}_k}$$
4. Re-normalize to prevent numerical underflow/overflow: $\mathbf{x}_{k+1} = \frac{A\mathbf{x}_k}{\|A\mathbf{x}_k\|}$.
5. Convergence is achieved when $|\lambda^{(k)} - \lambda^{(k-1)}| < \varepsilon$. The rate of convergence is governed by the ratio $\left|\frac{\lambda_2}{\lambda_1}\right|^k$.

#### Matrix Deflation (Wielandt Deflation)

To find subdominant eigenvalues without resolving the full eigenspectrum, deflate the dominant eigenmode from symmetric matrix $A$:

$$A' = A - \lambda_1 \frac{\mathbf{v}_1 \mathbf{v}_1^T}{\|\mathbf{v}_1\|^2}$$

- $A'$ has eigenvalues $0, \lambda_2, \lambda_3, \ldots, \lambda_n$, with the remaining eigenvectors preserved.
- Applying Power Iteration to $A'$ extracts $\lambda_2$ and $\mathbf{v}_2$.
- > ⚠ **Caution:** Numerical errors from the dominant mode accumulate into $A'$. Deflation beyond $\lambda_2$ is generally inadvisable without re-orthogonalization.

#### Tridiagonal Toeplitz Matrices

For an $n \times n$ tridiagonal matrix with diagonal $b$, superdiagonal $a$, and subdiagonal $c$:
$$\lambda_k = b + 2\sqrt{ac}\cos\left(\frac{k\pi}{n+1}\right), \qquad v_{i, k} = 2\left(\frac{c}{a}\right)^{\!i/2}\sin\left(\frac{i k \pi}{n+1}\right), \quad k, i \in \{1, \dots, n\}$$

- When $a = c = -1$ and $b = 2$ (the 1D discrete Laplacian stencil):
  The recurrence relation $(A\mathbf{v})_i = -v_{i-1} + 2v_i - v_{i+1} = \left[2 - 2\cos\left(\frac{k\pi}{n+1}\right)\right]v_i$ establishes that the formula with $+2\sqrt{ac}\cos(\dots)$ maps index $k$ to spatial mode $k' = n + 1 - k$. The dominant modes found by Power Iteration correspond to high-frequency spatial modes ($k' = 5, 4$) displaying alternating component signs across adjacent coordinates.

#### Code — `MyLib.py`

```python
def Power_Iteration(A, x0=None, tol=1e-6, max_iter=1000):
    """
    Finds dominant eigenvalue and normalized eigenvector of matrix A via Power Iteration.
    Uses Rayleigh quotient: lambda = (A*x . x) / (x . x).
    """
    n = len(A)
    x = [1.0] * n if x0 is None else [float(v) for v in x0]
    norm_x = math.sqrt(sum(v * v for v in x))
    x = [v / norm_x for v in x]
    
    lam_prev = 0.0
    for it in range(1, max_iter + 1):
        y = [sum(A[i][j] * x[j] for j in range(n)) for i in range(n)]
        lam = sum(y[i] * x[i] for i in range(n))  # Rayleigh quotient (x is unit vector)
        norm_y = math.sqrt(sum(v * v for v in y))
        if norm_y == 0:
            break
        x_new = [v / norm_y for v in y]
        if it > 1 and abs(lam - lam_prev) < tol:
            return lam, x_new, it
        lam_prev, x = lam, x_new
    return lam, x, max_iter

def Matrix_Deflation(A, lam, v):
    """Wielandt deflation for symmetric matrix: A' = A - lam * (v * v^T) / (v^T * v)."""
    n = len(A)
    norm_sq = sum(vi * vi for vi in v)
    A_def = [row[:] for row in A]
    for i in range(n):
        for j in range(n):
            A_def[i][j] -= lam * (v[i] * v[j]) / norm_sq
    return A_def
```

---

### 9 · Interpolation & Data Fitting

> 📄 [`CS-341-slides/datafit.pdf`](CS-341-slides/datafit.pdf)

Addresses reconstructing continuous functions from discrete measurements $f(x_i) = y_i$ (interpolation) and finding optimal mathematical models describing noisy experimental data (least-squares curve fitting).

#### Interpolation vs. Extrapolation

- **Interpolation:** Estimating $y(x)$ for $x$ strictly bounded within the data domain $[x_{\min}, x_{\max}]$.
- **Extrapolation:** Estimating $y(x)$ outside the domain $[x_{\min}, x_{\max}]$. Highly prone to massive errors and divergence unless backed by exact physical laws.

#### Lagrange Interpolating Polynomial

Given $N+1$ distinct points $(x_0, y_0), (x_1, y_1), \ldots, (x_N, y_N)$, there exists a unique polynomial $P_N(x)$ of degree at most $N$ passing through all points:

$$P_N(x) = \sum_{i=0}^N y_i L_i(x), \qquad L_i(x) = \prod_{\substack{k=0 \\ k \ne i}}^N \frac{x - x_k}{x_i - x_k}$$

- **Lagrange basis property:** $L_i(x_j) = \delta_{ij}$, ensuring $P_N(x_i) = y_i$ exactly.
- **Degree selection:** An order-$m$ polynomial requires $m+1$ points. If $m < N$, choosing the $m+1$ points closest to the evaluation point $x$ minimizes interpolation error.

#### Linear Regression & Least Squares

Given $N$ data points $(x_i, y_i)$ with measurement uncertainties $\sigma_i$ (default $\sigma_i = 1$), fit a straight line $f(x) = a_1 + a_2 x$ by minimizing the weighted sum of squared residuals:

$$\chi^2 = \sum_{i=1}^N \left(\frac{y_i - a_1 - a_2 x_i}{\sigma_i}\right)^{\!2}$$

Setting $\frac{\partial \chi^2}{\partial a_1} = 0$ and $\frac{\partial \chi^2}{\partial a_2} = 0$ yields the normal equations:

$$S = \sum_{i=1}^N \frac{1}{\sigma_i^2}, \quad S_x = \sum_{i=1}^N \frac{x_i}{\sigma_i^2}, \quad S_y = \sum_{i=1}^N \frac{y_i}{\sigma_i^2}, \quad S_{xx} = \sum_{i=1}^N \frac{x_i^2}{\sigma_i^2}, \quad S_{xy} = \sum_{i=1}^N \frac{x_i y_i}{\sigma_i^2}, \quad S_{yy} = \sum_{i=1}^N \frac{y_i^2}{\sigma_i^2}$$

$$\Delta = S S_{xx} - S_x^2$$

$$a_1 = \frac{S_{xx} S_y - S_x S_{xy}}{\Delta} \quad (\text{intercept}), \qquad a_2 = \frac{S S_{xy} - S_x S_y}{\Delta} \quad (\text{slope})$$

Parameter error estimates:
$$\sigma_{a_1} = \sqrt{\frac{S_{xx}}{\Delta}}, \qquad \sigma_{a_2} = \sqrt{\frac{S}{\Delta}}$$

#### Pearson's Correlation Coefficient ($r^2$)

Quantifies the linear association between $x$ and $y$:

$$r^2 = \frac{(S S_{xy} - S_x S_y)^2}{(S S_{xx} - S_x^2)(S S_{yy} - S_y^2)}, \qquad 0 \le r^2 \le 1$$

- $r^2 \to 1$: Strong linear correlation and superior fit quality.
- $r^2 \to 0$: Poor or nonexistent linear relationship.

#### Non-Linear Model Linearization

Many non-linear physical laws transform into linear forms via logarithms:

| Model | Formula | Linear Transformation | Linear Variables | Recovered Parameters |
|---|---|---|---|---|
| **Exponential** | $y = a e^{bx}$ | $\ln(y) = \ln(a) + bx$ | $Y = \ln(y), \ X = x$ | $a = e^{A_1}, \ b = A_2$ |
| **Exponential Decay** | $y = a e^{-bx}$ | $\ln(y) = \ln(a) - bx$ | $Y = \ln(y), \ X = x$ | $a = e^{A_1}, \ b = -A_2$ |
| **Power Law** | $y = a x^b$ | $\ln(y) = \ln(a) + b\ln(x)$ | $Y = \ln(y), \ X = \ln(x)$ | $a = e^{A_1}, \ b = A_2$ |
| **Logarithmic** | $y = a + b\ln(x)$ | $y = a + b\ln(x)$ | $Y = y, \ X = \ln(x)$ | $a = A_1, \ b = A_2$ |

The model producing a higher Pearson $r^2$ provides the statistically superior representation of the physical process.

#### Least Squares Polynomial Fitting & Normal Equations

Given $N$ data points $(x_i, y_i)$ with weights $w_i = 1/\sigma_i^2$ (default $w_i = 1$), fit a standard polynomial of degree $n$:

$$y = f(x) = \sum_{k=0}^n a_k x^k = a_0 + a_1 x + a_2 x^2 + \dots + a_n x^n$$

Minimizing the sum of squared residuals $R^2 = \sum_{i=1}^N w_i \left(y_i - \sum_{k=0}^n a_k x_i^k\right)^2$ with respect to each parameter $a_m$ ($\frac{\partial R^2}{\partial a_m} = 0$) produces the normal equations:

$$\sum_{k=0}^n a_k \left(\sum_{i=1}^N w_i x_i^{m+k}\right) = \sum_{i=1}^N w_i y_i x_i^m, \qquad m = 0, 1, \dots, n$$

In matrix form:

$$X \mathbf{a} = \mathbf{Y} \implies \mathbf{a} = X^{-1} \mathbf{Y}$$

where $X$ is an $(n+1) \times (n+1)$ **Hankel / normal matrix** and $\mathbf{Y}$ is the RHS vector:

$$X = \begin{pmatrix} \sum w_i & \sum w_i x_i & \dots & \sum w_i x_i^n \\ \sum w_i x_i & \sum w_i x_i^2 & \dots & \sum w_i x_i^{n+1} \\ \vdots & \vdots & \ddots & \vdots \\ \sum w_i x_i^n & \sum w_i x_i^{n+1} & \dots & \sum w_i x_i^{2n} \end{pmatrix}, \qquad \mathbf{Y} = \begin{pmatrix} \sum w_i y_i \\ \sum w_i y_i x_i \\ \vdots \\ \sum w_i y_i x_i^n \end{pmatrix}$$

#### Matrix Inversion: Cholesky vs. LU Decomposition

- **Symmetric Matrix:** $X_{mk} = X_{km}$, so $X$ is symmetric. For distinct data points, $X = \Phi^T W \Phi$ is symmetric positive definite.
- **Inversion via Cholesky Decomposition:** Solves $X \mathbf{x}_j = \mathbf{e}_j$ column-by-column using `Cholesky(A)` (storing the lower triangular factor $L$ with the `"CHOLESKY"` sentinel string in the upper triangle), packed via `Cholesky_compress(A)`, forward-solved with `forward()`, and back-substituted with $L^T = \text{Cholesky\_read(flat, upper=1)}$:
  $$L \mathbf{y} = \mathbf{e}_j \implies L^T \mathbf{x}_j = \mathbf{y}$$
- **Fallback to LU:** If the matrix is non-symmetric or fails positive-definiteness due to floating-point roundoff, it falls back to LU decomposition with partial pivoting (`LU(A)`).

#### Matrix Condition Number & Monomial Ill-Conditioning

The condition number $\kappa(X) = \|X\| \cdot \|X^{-1}\|$ measures sensitivity to roundoff errors and perturbation:
- $\kappa(X) \le 100$: Well-conditioned.
- $\kappa(X) \ge 10^4$: Ill-conditioned.

For standard polynomial bases $\{1, x, x^2, \dots, x^n\}$, entries grow rapidly as $x^{2n}$, and adjacent rows become nearly linearly dependent. The condition number explodes exponentially ($\kappa \sim 10^4$ for $n=3$, reaching $\sim 10^8 - 10^{16}$ for higher orders), exhausting 64-bit IEEE floating-point precision ($\sim 10^{-16}$) and destabilizing solutions.

#### Modified Chebyshev Polynomial Basis via 3-Term Recurrence

To circumvent ill-conditioning, orthogonal polynomial basis functions are used. For data scattered on $[0, 1]$, **shifted (modified) Chebyshev polynomials** $T_k^*(x) = T_k(2x - 1)$ are employed:

$$y = f(x) = \sum_{k=0}^n a_k \phi_k(x) = a_0 \phi_0(x) + a_1 \phi_1(x) + \dots + a_n \phi_n(x)$$

Explicit low-order polynomials:
$$\phi_0(x) = 1, \qquad \phi_1(x) = 2x - 1, \qquad \phi_2(x) = 8x^2 - 8x + 1, \qquad \phi_3(x) = 32x^3 - 48x^2 + 18x - 1$$

##### 3-Term Recurrence Relation
Generating functions require high-order derivatives, which are numerically unsuitable. Instead, $\phi_k(x)$ is computed efficiently via the 3-term recurrence relation:

$$\phi_0(x) = 1$$
$$\phi_1(x) = 2x - 1$$
$$\phi_{m+1}(x) = 2(2x - 1)\phi_m(x) - \phi_{m-1}(x) = (4x - 2)\phi_m(x) - \phi_{m-1}(x), \qquad m \ge 1$$

##### Why Chebyshev Basis Improves Conditioning
The modified Chebyshev normal matrix has entries:

$$X_{mk} = \sum_{i=1}^N \phi_m(x_i) \phi_k(x_i)$$

Due to orthogonality over the domain, cross-terms average out ($\sum \phi_m \phi_k \approx \delta_{mk}$), making $X$ strongly diagonally dominant.

**Experimental comparison on `assign9fit.txt` ($n = 3$):**
- Standard Polynomial Condition Number: $\kappa(X_{\text{poly}}) \approx 1.2105 \times 10^4$
- Modified Chebyshev Condition Number: $\kappa(X_{\text{cheb}}) \approx 3.8561$
- **Improvement:** $> 3139\times$ reduction in condition number, providing exceptional numerical stability for Cholesky inversion while achieving identical goodness of fit ($r^2 = 0.9724$).

#### Code — `MyLib.py`

```python
def is_symmetric_matrix(A, tol=1e-9):
    """Checks whether square matrix A is symmetric within given tolerance."""
    n = len(A)
    for i in range(n):
        for j in range(i + 1, n):
            if not math.isclose(A[i][j], A[j][i], rel_tol=tol, abs_tol=tol):
                return False
    return True

def Solve_Linear(A, B):
    """Solves A*x = B via Cholesky (if symmetric positive definite) or LU with pivoting."""
    n = len(A)
    A_copy = [row[:] for row in A]
    if is_symmetric_matrix(A_copy):
        try:
            Cholesky(A_copy)
            flat = Cholesky_compress(A_copy)
            y = forward(flat, B)
            return backward(Cholesky_read(flat, upper=1), y)
        except Exception:
            pass
    A_lu, perm = LU(A_copy)
    y = forward(A_lu, B, perm=perm)
    return backward(A_lu, y)

def Invert_Matrix(A):
    """Inverts matrix A column-by-column using Cholesky (if symmetric) or LU decomposition."""
    n = len(A)
    inv = [[0.0] * n for _ in range(n)]
    if is_symmetric_matrix(A):
        try:
            A_copy = [row[:] for row in A]
            Cholesky(A_copy)
            flat = Cholesky_compress(A_copy)
            LT = Cholesky_read(flat, upper=1)
            for j in range(n):
                e = [1.0 if i == j else 0.0 for i in range(n)]
                inv_col = backward(LT, forward(flat, e))
                for i in range(n):
                    inv[i][j] = inv_col[i]
            return inv
        except Exception:
            pass
    A_lu, perm = LU([row[:] for row in A])
    for j in range(n):
        e = [1.0 if i == j else 0.0 for i in range(n)]
        inv_col = backward(A_lu, forward(A_lu, e, perm=perm))
        for i in range(n):
            inv[i][j] = inv_col[i]
    return inv

def Modified_Chebyshev_Poly(k, x):
    """Evaluates modified Chebyshev basis polynomial phi_k(x) via 3-term recurrence."""
    if k == 0: return 1.0
    if k == 1: return 2.0 * x - 1.0
    p_prev2, p_prev1 = 1.0, 2.0 * x - 1.0
    factor = 4.0 * x - 2.0
    for _ in range(2, k + 1):
        p_curr = factor * p_prev1 - p_prev2
        p_prev2, p_prev1 = p_prev1, p_curr
    return p_prev1

def Poly_Fit(dataset, order, sigmas=None):
    """Fits standard polynomial y = sum_{k=0}^n a_k * x^k via normal equations and Cholesky/LU."""
    pts = [[float(row[0]), float(row[1])] for row in dataset]
    n_pts, m_dim = len(pts), order + 1
    w = [1.0] * n_pts if sigmas is None else [1.0 / (s * s) for s in sigmas]
    X = [[sum(w[i] * (pts[i][0] ** (m + k)) for i in range(n_pts)) for k in range(m_dim)] for m in range(m_dim)]
    Y = [sum(w[i] * pts[i][1] * (pts[i][0] ** m) for i in range(n_pts)) for m in range(m_dim)]
    coeffs = Solve_Linear(X, Y)
    X_inv = Invert_Matrix(X)
    cond_num = Matrix_Condition_Number(X)
    y_mean = sum(p[1] for p in pts) / n_pts
    ss_tot = sum((p[1] - y_mean) ** 2 for p in pts)
    ss_res = sum((p[1] - sum(coeffs[k] * (p[0] ** k) for k in range(m_dim))) ** 2 for p in pts)
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0
    return coeffs, X, Y, X_inv, cond_num, r2

def Chebyshev_Fit(dataset, order, sigmas=None):
    """Fits data using modified Chebyshev basis functions via 3-term recurrence and Cholesky/LU."""
    pts = [[float(row[0]), float(row[1])] for row in dataset]
    n_pts, m_dim = len(pts), order + 1
    w = [1.0] * n_pts if sigmas is None else [1.0 / (s * s) for s in sigmas]
    X = [[sum(w[i] * Modified_Chebyshev_Poly(m, pts[i][0]) * Modified_Chebyshev_Poly(k, pts[i][0]) for i in range(n_pts)) for k in range(m_dim)] for m in range(m_dim)]
    Y = [sum(w[i] * pts[i][1] * Modified_Chebyshev_Poly(m, pts[i][0]) for i in range(n_pts)) for m in range(m_dim)]
    coeffs = Solve_Linear(X, Y)
    X_inv = Invert_Matrix(X)
    cond_num = Matrix_Condition_Number(X)
    y_mean = sum(p[1] for p in pts) / n_pts
    ss_tot = sum((p[1] - y_mean) ** 2 for p in pts)
    ss_res = sum((p[1] - sum(coeffs[k] * Modified_Chebyshev_Poly(k, p[0]) for k in range(m_dim))) ** 2 for p in pts)
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0
    return coeffs, X, Y, X_inv, cond_num, r2
```

```python
def Lagrange_Interp(dataset, x_eval, order=None):
    """Lagrange interpolation of order m (uses m+1 nearest points if order specified)."""
    pts = [[float(row[0]), float(row[1])] for row in dataset]
    if order is not None and order + 1 < len(pts):
        pts.sort(key=lambda p: abs(p[0] - x_eval))
        pts = pts[:order + 1]
    m = len(pts)
    y_eval = 0.0
    for i in range(m):
        term = pts[i][1]
        for k in range(m):
            if k != i:
                term *= (x_eval - pts[k][0]) / (pts[i][0] - pts[k][0])
        y_eval += term
    return y_eval

def Linear_Regression(x, y, sigmas=None):
    """Weighted/unweighted linear fit y = a1 + a2*x, returning a1, a2, Pearson r^2, sig_a1, sig_a2."""
    n = len(x)
    w = [1.0] * n if sigmas is None else [1.0 / (s * s) for s in sigmas]
    S = sum(w)
    Sx, Sy = sum(w[i] * x[i] for i in range(n)), sum(w[i] * y[i] for i in range(n))
    Sxx = sum(w[i] * x[i] * x[i] for i in range(n))
    Sxy = sum(w[i] * x[i] * y[i] for i in range(n))
    Syy = sum(w[i] * y[i] * y[i] for i in range(n))
    Delta = S * Sxx - Sx * Sx
    a1 = (Sxx * Sy - Sx * Sxy) / Delta
    a2 = (S * Sxy - Sx * Sy) / Delta
    sig_a1, sig_a2 = math.sqrt(Sxx / Delta), math.sqrt(S / Delta)
    r2 = ((S * Sxy - Sx * Sy) ** 2) / ((S * Sxx - Sx * Sx) * (S * Syy - Sy * Sy))
    return a1, a2, r2, sig_a1, sig_a2
```

---

## 📚 Reference Book

**Alex Gezerlis — *Numerical Methods in Physics with Python***

A superb companion to this course — covers all the methods above with rigorous derivations and Python implementations.

---

<div align="center">

*Notes compiled from lecture slides by Prof. Subhasish Basak, PHY341/745, NISER.*

</div>
