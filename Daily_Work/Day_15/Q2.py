import math
import matplotlib.pyplot as plt
from MyLib import *

L = 2.0
T = 4.0
Nx = 20
Nt = 2000

hx = L / Nx
ht = T / Nt

alpha = ht / (hx ** 2)

print("1D Heat Equation - Explicit Scheme")
print("Length of rod L =", L)
print("Total time T =", T)
print("Number of spatial intervals Nx =", Nx)
print("Number of time steps Nt =", Nt)
print("Spatial step size hx =", hx)
print("Time step size ht =", ht)
print("Stability ratio alpha = ht / (hx^2) =", alpha)

if alpha <= 0.5:
	print("Stability condition (alpha <= 0.5) is satisfied.")
else:
	print("Warning: Stability condition (alpha <= 0.5) is violated.")
print()

# Spatial grid points
xs = []
for i in range(Nx + 1):
	position = i * hx
	xs.append(position)

# Initial condition at t = 0: u(x, 0) = 20 * |sin(pi * x)|
u = []
for i in range(Nx + 1):
	x_val = xs[i]
	temp_val = 20.0 * abs(math.sin(math.pi * x_val))
	u.append(temp_val)

u[0] = 0.0
u[Nx] = 0.0

target_steps = [0, 10, 20, 50, 100, 200, 500, 1000]

saved_profiles = {}
saved_profiles[0] = list(u)

# Time evolution loop using explicit scheme
max_step = max(target_steps)

for step in range(1, max_step + 1):
	u_new = [0.0] * (Nx + 1)
	
	for i in range(1, Nx):
		u_new[i] = alpha * (u[i + 1] + u[i - 1]) + (1.0 - 2.0 * alpha) * u[i]
	
	# Boundary conditions at x = 0 and x = L
	u_new[0] = 0.0
	u_new[Nx] = 0.0
	
	u = u_new
	
	if step in target_steps:
		saved_profiles[step] = list(u)

print(f"{'Step':<8} {'Time t':<10} {'Max Temp (C)':<15} {'Center Temp (x=1.0)':<20}")
for step in target_steps:
	current_time = step * ht
	profile = saved_profiles[step]
	max_temp = max(profile)
	center_temp = profile[Nx // 2]
	print(f"{step:<8} {current_time:<10.3f} {max_temp:<15.4f} {center_temp:<20.4f}")

plt.figure()

for step in target_steps:
	current_time = step * ht
	plt.plot(xs, saved_profiles[step], marker="o", markersize=3, label=f"step {step} (t = {current_time:.3f})")

plt.xlabel("Position along the rod x")
plt.ylabel("Temperature u(x, t) in C")
plt.title("Q2: Temperature Profile across the Rod (Explicit Scheme)")
plt.legend()
plt.grid(True)
plt.show()
