from MyLib import *
import matplotlib.pyplot as plt

data = []
with open("assign9fit.txt", "r") as file:
	for line in file:
		line = line.strip()
		if not line or line.startswith("#"):
			continue
		parts = line.split()
		if len(parts) >= 2:
			data.append([float(parts[0]), float(parts[1])])

print("Question 3: Cubic Least Squares Polynomial Fit")
print(f"Loaded {len(data)} data points from assign9fit.txt\n")

order = 3
coeffs, X, Y, X_inv, cond_num, r2 = Poly_Fit(data, order=order)

print("Normal Matrix X (Hankel Form):")
for row in X:
	print("  " + "  ".join(f"{val:14.6e}" for val in row))
print()

print("RHS Vector Y:")
for val in Y:
	print(f"  {val:14.6e}")
print()

print("Inverse Matrix X^-1 (via Cholesky Decomposition):")
for row in X_inv:
	print("  " + "  ".join(f"{val:14.6e}" for val in row))
print()

print("Fitted Coefficients for y = a0 + a1*x + a2*x^2 + a3*x^3:")
for k, ak in enumerate(coeffs):
	print(f"  a{k} = {ak:14.6f}")
print()

print(f"Condition Number kappa(X): {cond_num:.6e}")
print(f"Coefficient of Determination (Pearson r^2): {r2:.6f}")

# Plotting experimental data and the cubic fit
xs = [p[0] for p in data]
ys = [p[1] for p in data]

x_min, x_max = min(xs), max(xs)
span = x_max - x_min
x_dense = [x_min + i * (span / 200.0) for i in range(201)]
y_fit = [sum(coeffs[k] * (x ** k) for k in range(order + 1)) for x in x_dense]

plt.figure(figsize=(10, 6))
plt.scatter(xs, ys, color="black", zorder=5, label="Experimental Data (assign9fit.txt)")
plt.plot(x_dense, y_fit, "r-", linewidth=2, label=f"Cubic Polynomial Fit ($r^2={r2:.4f}$)")
plt.xlabel("x", fontsize=12)
plt.ylabel("y", fontsize=12)
plt.title("Question 3: Cubic Least Squares Polynomial Fit", fontsize=14)
plt.legend(fontsize=11)
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.show()
