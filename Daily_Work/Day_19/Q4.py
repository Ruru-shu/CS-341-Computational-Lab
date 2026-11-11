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

print("Question 4: Modified Chebyshev Polynomial Least Squares Fit")
print(f"Loaded {len(data)} data points from assign9fit.txt\n")

order = 3

# Modified Chebyshev Fit
c_cheb, X_cheb, Y_cheb, X_inv_cheb, cond_cheb, r2_cheb = Chebyshev_Fit(data, order=order)

# Standard Polynomial Fit for comparison
c_poly, X_poly, Y_poly, X_inv_poly, cond_poly, r2_poly = Poly_Fit(data, order=order)

print("Modified Chebyshev Normal Matrix X_cheb:")
for row in X_cheb:
	print("  " + "  ".join(f"{val:14.6e}" for val in row))
print()

print("RHS Vector Y_cheb:")
for val in Y_cheb:
	print(f"  {val:14.6e}")
print()

print("Inverse Matrix X_cheb^-1 (via Cholesky Decomposition):")
for row in X_inv_cheb:
	print("  " + "  ".join(f"{val:14.6e}" for val in row))
print()

print("Fitted Coefficients for y = a0*phi0 + a1*phi1 + a2*phi2 + a3*phi3:")
for k, ak in enumerate(c_cheb):
	print(f"  a{k} = {ak:14.6f}")
print()

print(f"Condition Number kappa(X_cheb): {cond_cheb:.6e}")
print(f"Coefficient of Determination (Pearson r^2): {r2_cheb:.6f}")
print("=" * 70)
print("COMPARISON & RELATIVE MERIT OF BASIS FUNCTIONS:")
print("=" * 70)
print(f"Standard Polynomial Condition Number kappa(X_poly) : {cond_poly:.6e}")
print(f"Modified Chebyshev Condition Number kappa(X_cheb)  : {cond_cheb:.6e}")
ratio = cond_poly / cond_cheb if cond_cheb > 0 else float("inf")
print(f"Condition Number Ratio (kappa_poly / kappa_cheb)  : {ratio:.2f}x")
print()
print("Explanation:")
print("1. Standard polynomial monomials {1, x, x^2, x^3} lead to a Hankel normal matrix")
print("   where elements grow rapidly (sum x^(m+k)), making rows nearly linearly dependent.")
print("   This causes severe ill-conditioning and makes matrix inversion susceptible to floating-point roundoff.")
print("2. In contrast, modified Chebyshev polynomials {phi_0, phi_1, phi_2, phi_3} form an")
print("   orthogonal basis on [0, 1]. The inner products sum phi_m(x_i)*phi_k(x_i) approx delta_mk,")
print("   yielding a well-conditioned, diagonally dominant matrix.")
print("3. Consequently, the condition number for Chebyshev basis is orders of magnitude smaller,")
print("   ensuring robust, numerically stable matrix inversion via Cholesky decomposition.")
print("=" * 70)

# Plotting data and both fits
xs = [p[0] for p in data]
ys = [p[1] for p in data]

x_min, x_max = min(xs), max(xs)
span = x_max - x_min
x_dense = [x_min + i * (span / 200.0) for i in range(201)]

y_fit_poly = [sum(c_poly[k] * (x ** k) for k in range(order + 1)) for x in x_dense]
y_fit_cheb = [sum(c_cheb[k] * Modified_Chebyshev_Poly(k, x) for k in range(order + 1)) for x in x_dense]

plt.figure(figsize=(10, 6))
plt.scatter(xs, ys, color="black", zorder=5, label="Experimental Data (assign9fit.txt)")
plt.plot(x_dense, y_fit_poly, "r--", linewidth=2, label=f"Standard Cubic Fit ($\\kappa={cond_poly:.1e}$)")
plt.plot(x_dense, y_fit_cheb, "b-", alpha=0.7, linewidth=2, label=f"Modified Chebyshev Fit ($\\kappa={cond_cheb:.1e}$)")
plt.xlabel("x", fontsize=12)
plt.ylabel("y", fontsize=12)
plt.title("Least Squares Polynomial vs. Modified Chebyshev Fitting", fontsize=14)
plt.legend(fontsize=11)
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.show()
