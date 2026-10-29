import math
import matplotlib.pyplot as plt
from MyLib import *

# Problem parameters from Question 2:
# ut = a * uxx + b * ux + c * u
# a = 2, b = -1, c = 1, on 0 <= x <= 1
a = 2.0
b = -1.0
c = 1.0

L = 1.0
T = 1.0

Nx = 50
Nt = 1000

# Initial condition: u(x, 0) = 2 + x^2
def initial_condition(x):
	return 2.0 + x * x

# Boundary conditions: u(0, t) = 2, u(1, t) = 3
u0 = 2.0
uL = 3.0

print("1D Advection-Diffusion-Reaction PDE - Crank-Nicolson Scheme")
print("Domain: 0 <= x <= 1, 0 <= t <= 1")
print(f"Parameters: a = {a}, b = {b}, c = {c}")
print(f"Boundary conditions: u(0, t) = {u0}, u(1, t) = {uL}")
print(f"Number of spatial intervals Nx = {Nx}")
print(f"Number of time steps Nt = {Nt}")

xs, ts, u, A, B, alpha, beta, gamma = Crank_Nicolson_ADR(
	a=a, b=b, c=c, L=L, T=T, Nx=Nx, Nt=Nt,
	f_init=initial_condition, u0=u0, uL=uL
)

hx = L / Nx
ht = T / Nt

print(f"Spatial step size hx = {hx:.4f}")
print(f"Time step size ht = {ht:.6f}")
print(f"Parameters: alpha = {alpha:.6f}, beta = {beta:.6f}, gamma = {gamma:.6f}")
print()

# Sample time steps to display
sample_times = [0.0, 0.05, 0.1, 0.25, 0.5, 1.0]
sample_steps = [int(round(st / ht)) for st in sample_times]

print(f"{'Time t':<10} {'u(0.0)':<10} {'u(0.25)':<12} {'u(0.5)':<12} {'u(0.75)':<12} {'u(1.0)':<10}")
for st, step in zip(sample_times, sample_steps):
	idx_25 = int(round(0.25 * Nx))
	idx_50 = int(round(0.50 * Nx))
	idx_75 = int(round(0.75 * Nx))
	print(f"{st:<10.3f} {u[step][0]:<10.4f} {u[step][idx_25]:<12.4f} {u[step][idx_50]:<12.4f} {u[step][idx_75]:<12.4f} {u[step][Nx]:<10.4f}")

# Figure 1: Contour Plot of u(x, t)
plt.figure()
contour_fill = plt.contourf(xs, ts, u, levels=40, cmap="viridis")
cbar = plt.colorbar(contour_fill)
cbar.set_label("Solution u(x, t)")
plt.contour(xs, ts, u, levels=12, colors="black", linewidths=0.5, alpha=0.6)
plt.xlabel("Position x")
plt.ylabel("Time t")
plt.title("Q2: Contour Plot of u(x, t) (Crank-Nicolson Scheme)")
plt.grid(True, linestyle="--", alpha=0.5)

# Figure 2: Temperature / Solution profiles across rod at different times
plt.figure()
for st, step in zip(sample_times, sample_steps):
	plt.plot(xs, u[step], marker="o", markersize=2, label=f"t = {st:.2f}")

plt.xlabel("Position x")
plt.ylabel("u(x, t)")
plt.title("Q2: Solution Profiles across Rod (Crank-Nicolson Scheme)")
plt.legend()
plt.grid(True)

plt.show()
