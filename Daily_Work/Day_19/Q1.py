from MyLib import *

data = []
with open("assign9fit.txt", "r") as file:
	for line in file:
		line = line.strip()
		if not line or line.startswith("#"):
			continue
		parts = line.split()
		if len(parts) >= 2:
			data.append([float(parts[0]), float(parts[1])])

print("Question 1: Library Function for Standard Polynomial Least Square Fitting")
print(f"Total data points loaded from assign9fit.txt: {len(data)}")
print("-" * 65)

for order in [1, 2, 3]:
	coeffs, X, Y, X_inv, cond_num, r2 = Poly_Fit(data, order=order)
	print(f"Order n = {order}:")
	print(f"  Coefficients [a0, ..., a{order}]: {[round(c, 6) for c in coeffs]}")
	print(f"  Condition Number kappa(X): {cond_num:.6e}")
	print(f"  Coefficient of determination r^2: {r2:.6f}")
	print("-" * 65)
