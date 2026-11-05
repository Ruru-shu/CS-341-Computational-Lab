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
│   └── eigv.pdf            # Eigenvalues & Eigenvectors (Power Iteration & Deflation)
│
├── MyLib.py                # Core library: PRNG, MyComplex
├── MyLib_Advanced.py       # Extended library: Solve, Newton-Raphson,
│                           #   Fixed-Point, Deriv, GaussPRNG
├── PRNG.py                 # Standalone PRNG module
├── PRNG_Matrices.py        # PRNG + matrix utilities
│
├── Daily_Work/
│   ├── Day_0/ … Day_17/    # Daily class work, scripts, and outputs
│   └── (Day_17/MyLib.py)   # Most complete library snapshot
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

## 📚 Reference Book

**Alex Gezerlis — *Numerical Methods in Physics with Python***

A superb companion to this course — covers all the methods above with rigorous derivations and Python implementations.

---

<div align="center">

*Notes compiled from lecture slides by Prof. Subhasish Basak, PHY341/745, NISER.*

</div>
