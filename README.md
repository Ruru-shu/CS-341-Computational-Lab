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
│   └── integ.pdf           # Numerical integration
│
├── MyLib.py                # Core library: PRNG, MyComplex
├── MyLib_Advanced.py       # Extended library: Solve, Newton-Raphson,
│                           #   Fixed-Point, Deriv, GaussPRNG
├── PRNG.py                 # Standalone PRNG module
├── PRNG_Matrices.py        # PRNG + matrix utilities
│
├── Day_0/ … Day_9/         # Daily class work and assignments
│   └── (Day_9/MyLib.py)    # Most complete library snapshot
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

Error bound: $\displaystyle\left|\mathcal{M}_N - \mathcal{I}\right| \le \frac{(b-a)^3}{24N^2}\,|f''|_{\max}$

#### Trapezoidal rule

Connect consecutive $f(x_n)$ with straight lines:

$$\mathcal{T}_N = \frac{h}{2}\!\left[f(x_0) + 2f(x_1) + \cdots + 2f(x_{N-1}) + f(x_N)\right]$$

Error bound: $\displaystyle\left|\mathcal{T}_N - \mathcal{I}\right| \le \frac{(b-a)^3}{12N^2}\,|f''|_{\max}$

> ⚠ Surprisingly, the trapezoidal rule can be *less* accurate than midpoint for strictly concave or convex functions.

#### Simpson's 1/3 rule

Uses parabolic interpolation over pairs of sub-intervals (**even** $N$ required):

$$\mathcal{S}_N = \frac{h}{3}\!\left[f(x_0) + 4f(x_1) + 2f(x_2) + 4f(x_3) + \cdots + 4f(x_{N-1}) + f(x_N)\right]$$

Weights follow the pattern $1, 4, 2, 4, 2, \ldots, 4, 1$.

Error bound: $\displaystyle\left|\mathcal{S}_N - \mathcal{I}\right| \le \frac{(b-a)^5}{180N^4}\,|f''''|_{\max}$

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

## 📚 Reference Book

**Alex Gezerlis — *Numerical Methods in Physics with Python***

A superb companion to this course — covers all the methods above with rigorous derivations and Python implementations.

---

<div align="center">

*Notes compiled from lecture slides by Prof. Subhasish Basak, PHY341/745, NISER.*

</div>
