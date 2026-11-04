from MyLib import *

A = []
with open("A2.txt", "r") as file:
	for line in file:
		A.append(list(map(float, line.strip().split())))

print("Matrix A:")
for row in A:
	print("  " + "  ".join(f"{val:8.4f}" for val in row))
print()

eigenvalue, eigenvector, iters = Power_Iteration(A, tol=1e-6)

print("Power Iteration Method Results:")
print(f"Dominant Eigenvalue: {eigenvalue:.6f}")
print("Corresponding Eigenvector (normalized):")
print("  v =", [round(x, 6) for x in eigenvector])
print(f"Iterations to converge (tol=1e-6): {iters}")
