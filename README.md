#2026-08-05

Instructor: Prof. Subhasis Basak

Name - Omm Biswal
Roll - 2411122
# Course overview

Enter class asap before the proff
Assignments basically kinda daily
All communications over Google Classroom
Midsem exam the first class after Midsem
Not >7 kids gonna get AA

Assignments - 35
Midsem - 15
DIY Project - 20
Endsem - 30

# Alex Gerzerliks - Numerical Methods in Physics with Python
---
## Algorithms to be covered

![[Pasted image 20260912192335.png]]

RNG
Gauss-Jordan Elimination
LU decomposition
Jacobi Gauss Seidel - iterative of the above
Roots finding via Regula Falsi, Netwon-Raphson, Laguerre (polynomials)
Integration - Midpoint, Trapezoidal, Simpson, Monte Carlo
ODE - Euler, Predictor-Corrector, Runge-Kutta (both initial and boundary value problems)
Least square fitting

If time did allow, PDE & Eigenvalue problem, ML regression and classification

---
#2026-08-12

# Random Numbers

Basic Goal: Writing own PRNG and use it for all assignments, exams, DIY etc.
$x_{n+1} = cx_{n}(1-x_{n})$

Basic test for randomness: Correlation, moments
Advanced test: chi-square, Kolmogorov-Smirnov

Ideally Random numbers generated have n correlations and error statistical i.e, scale as $1/\sqrt{N}$

Linear Congruential Generator (LCG)

Use for Radioactive Decay

Tmr about exponential distribution

---

#19-08-2026

# LU Decomposition

Lower-Upper Decomposition, s.t.

L . U = A for Lower Triangular x Upper Triangular Matrix = A, original matrix, for
A.x = B

//Dolittle LU factorization

// Read matrix from a file
Overwriting.. 

# Gauss_Seidel and Jacobi

Make a temp array in Gauss_Seidel to compare
You must use magnitude of the whole vector to check for tolerance

Successive Over-Relaxation method

![[Pasted image 20260827153913.png]]


# Finding Roots

## Bisection Method

Basically using bisection and bracketting via IVT (Intermediate Value Theorem)

Bracketting fails if we get f(b) and f(a) with the same sign. We can something about that, mm
> [!note]-
> If both +ve or -ve, we can shift one of them by a value $\beta(b-a)$; $\beta=[0,1]$ depending on choice

If we get opposite sign, and f(x) in continous in that interval \[a,b\] then we keep using Binary search to lower the bracketting interval until b-a has a tolerance level of e-6 or something. And then we declare the ans to be avg(a, b).

## Regula-Falsi (False Position)

Interpolation and convergence faster than Bisection.
Bracketting but once we [""|bracket it], we can find slope of the function and find the slope of the line joining f(a) and f(b), and find the point that crosses the x-axis.
Then, you get that c/a c.
Then you replace a or b with this c and repeat.
Then after you get that tolerance is lower than a certain value, you say $c_n$ is the root.


## Laguerre
### Midpoint

### Trapezoidal
