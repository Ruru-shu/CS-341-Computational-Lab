import math
import matplotlib.pyplot as plt
from MyLib import *

# Read data table from data_q4.txt using open()
data = []
with open("data_q4.txt", "r") as file:
	for line in file:
		parts = line.strip().split()
		if parts:
			data.append([float(parts[0]), float(parts[1])])

xs = [row[0] for row in data]
ys = [row[1] for row in data]
n = len(xs)

print("Question 4: Non-linear Model Fitting via Linear Regression")
print("Data Points:")
print(f"{'x':>8} {'y':>8}")
for x_val, y_val in zip(xs, ys):
	print(f"{x_val:8.2f} {y_val:8.2f}")
print()

# Model 1: Power Law y = a * x^b
# Linear form: ln(y) = ln(a) + b * ln(x) -> Y = A1 + A2 * X
X_pow = [math.log(x) for x in xs]
Y_pow = [math.log(y) for y in ys]
A1_pow, A2_pow, r2_pow, sig_A1_pow, sig_A2_pow = Linear_Regression(X_pow, Y_pow)
a_pow = math.exp(A1_pow)
b_pow = A2_pow

print("1. Power Law Model: y = a * x^b")
print(f"   Transformed Equation: ln(y) = {A1_pow:.4f} + ({A2_pow:.4f}) * ln(x)")
print(f"   Parameters: a = {a_pow:.4f}, b = {b_pow:.4f}")
print(f"   Fitted Model: y = {a_pow:.4f} * x^({b_pow:.4f})")
print(f"   Pearson's r^2 = {r2_pow:.6f}\n")

# Model 2: Exponential y = a * e^(-b * x)
# Linear form: ln(y) = ln(a) - b * x -> Y = A1 + A2 * X where A2 = -b
X_exp = xs[:]
Y_exp = [math.log(y) for y in ys]
A1_exp, A2_exp, r2_exp, sig_A1_exp, sig_A2_exp = Linear_Regression(X_exp, Y_exp)
a_exp = math.exp(A1_exp)
b_exp = -A2_exp

print("2. Exponential Model: y = a * e^(-b * x)")
print(f"   Transformed Equation: ln(y) = {A1_exp:.4f} + ({A2_exp:.4f}) * x")
print(f"   Parameters: a = {a_exp:.4f}, b = {b_exp:.4f}")
print(f"   Fitted Model: y = {a_exp:.4f} * e^(-{b_exp:.4f} * x)")
print(f"   Pearson's r^2 = {r2_exp:.6f}\n")

# Model Comparison
print("=" * 60)
print(f"{'Model':<20} {'Fitted Equation':<28} {'Pearson r^2':<12}")
print("-" * 60)
print(f"{'Power Law':<20} {f'y = {a_pow:.2f} * x^({b_pow:.2f})':<28} {r2_pow:<12.6f}")
print(f"{'Exponential':<20} {f'y = {a_exp:.2f} * e^(-{b_exp:.3f}x)':<28} {r2_exp:<12.6f}")
print("=" * 60)

if r2_pow > r2_exp:
	print(f"Conclusion: The Power Law model provides a better fit (r^2 = {r2_pow:.6f} > {r2_exp:.6f}).")
else:
	print(f"Conclusion: The Exponential model provides a better fit (r^2 = {r2_exp:.6f} > {r2_pow:.6f}).")

# Generate smooth curves for plotting
x_dense = [xs[0] + i * (xs[-1] - xs[0]) / 200.0 for i in range(201)]
y_pow_dense = [a_pow * (x ** b_pow) for x in x_dense]
y_exp_dense = [a_exp * math.exp(-b_exp * x) for x in x_dense]

# Plot
plt.figure(figsize=(9, 6))
plt.scatter(xs, ys, color="black", s=45, zorder=5, label="Data Points")
plt.plot(x_dense, y_pow_dense, color="crimson", linewidth=2.0, label=f"Power Law: $y = {a_pow:.2f} x^{{{b_pow:.2f}}}$ ($r^2 = {r2_pow:.4f}$)")
plt.plot(x_dense, y_exp_dense, color="royalblue", linestyle="--", linewidth=2.0, label=f"Exponential: $y = {a_exp:.2f} e^{{-{b_exp:.3f}x}}$ ($r^2 = {r2_exp:.4f}$)")

plt.xlabel("x", fontsize=12)
plt.ylabel("y", fontsize=12)
plt.title("Q4: Comparison of Power Law vs Exponential Fits", fontsize=13)
plt.legend(fontsize=11)
plt.grid(True, linestyle="--", alpha=0.6)

plt.show()
