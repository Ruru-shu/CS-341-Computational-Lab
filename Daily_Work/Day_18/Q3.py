from MyLib import *

data = []
with open("data_q3.txt", "r") as file:
	for line in file:
		parts = line.strip().split()
		if parts:
			data.append([float(parts[0]), float(parts[1])])

print("Input Data Table:")
print(f"{'x':>8} {'y(x)':>10}")
for x_val, y_val in data:
	print(f"{x_val:8.2f} {y_val:10.2f}")
print()

x_target = 6.7
n_points = len(data)

print(f"Target point for estimation: x = {x_target}")
print("-" * 50)
print(f"{'Order':>6} {'Degree':>8} {'Points Used':>14} {'Estimate y(6.7)':>18}")

for order in range(1, n_points):
	est = Lagrange_Interp(data, x_target, order=order)
	pts_used = order + 1
	print(f"{order:6d} {order:8d} {pts_used:14d} {est:18.6f}")

# Full polynomial interpolation (order = n - 1 = 4)
best_estimate = Lagrange_Interp(data, x_target, order=n_points - 1)
print("-" * 50)
print(f"Best Estimate for y(x = {x_target}): {best_estimate:.6f}")
