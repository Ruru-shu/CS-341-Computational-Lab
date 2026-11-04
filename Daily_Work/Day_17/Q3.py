from MyLib import *

A = []
with open("assign8.txt", "r") as file:
	for line in file:
		A.append(list(map(float, line.strip().split())))

n = len(A)
print(f"Matrix A ({n}x{n} tridiagonal)")
for row in A:
	print("  " + "  ".join(f"{val:4.1f}" for val in row))
print()

a = -1.0
c = -1.0
b = 2.0

print("Analytical values given by the formula:")
print("Formula: lambda_k = b + 2*sqrt(ac)*cos(k*pi/(n+1))")
print("         v_{i,k} = 2*(c/a)^(i/2)*sin(i*k*pi/(n+1))")
print(f"Parameters: a = {a}, c = {c}, b = {b}, n = {n}\n")

analytical_vals = {}
analytical_vecs = {}
analytical_vecs_norm = {}

for k in range(1, n + 1):
	lam_k = b + 2.0 * math.sqrt(a * c) * math.cos(k * math.pi / (n + 1))
	# Note: c/a = (-1)/(-1) = 1, so (c/a)^(i/2) = 1
	v_k = [2.0 * math.sin(i * k * math.pi / (n + 1)) for i in range(1, n + 1)]
	norm_v = math.sqrt(sum(x * x for x in v_k))
	v_k_norm = [x / norm_v for x in v_k]

	analytical_vals[k] = lam_k
	analytical_vecs[k] = v_k
	analytical_vecs_norm[k] = v_k_norm
	print(f"k = {k}: lambda_{k} = {lam_k:.6f}")
	print(f"       v_{k} (raw)        = {[round(x, 4) for x in v_k]}")
	print(f"       v_{k} (normalized) = {[round(x, 4) for x in v_k_norm]}")

print("\n" + "=" * 65 + "\n")

# Numerical determination using Power Iteration method
# 1st dominant eigenvalue & eigenvector:
lam1_num, v1_num, iters1 = Power_Iteration(A, x0=[1.0] * n, tol=1e-8)

# Deflate matrix A using dominant eigenpair:
# A' = A - lambda_1 * (v1 v1^T)
A_def = Matrix_Deflation(A, lam1_num, v1_num)

# 2nd dominant eigenvalue & eigenvector:
# Initial guess chosen to have non-zero projection along the 2nd mode (e.g. anti-symmetric or monotonic)
lam2_num, v2_num, iters2 = Power_Iteration(A_def, x0=[1.0, 2.0, 3.0, 4.0, 5.0], tol=1e-8)

print("Power Iteration Method Numerical Results:")
print(f"1st Dominant Eigenvalue: {lam1_num:.6f}  (iters: {iters1})")
print(f"1st Eigenvector: {[round(x, 6) for x in v1_num]}")
print(f"2nd Dominant Eigenvalue: {lam2_num:.6f}  (iters: {iters2})")
print(f"2nd Eigenvector: {[round(x, 6) for x in v2_num]}")

print("\n" + "=" * 65 + "\n")

# Verification & Comparison
print("Verification with Analytical Expressions:")
print(f"Largest numerical eigenvalue:     {lam1_num:.6f}")
print(f"Analytical eigenvalue for k=1:    {analytical_vals[1]:.6f}")
print(f"Relative error:                   {abs(lam1_num - analytical_vals[1]) / analytical_vals[1]:.2e}")
print()
print(f"Second largest numerical eigenvalue: {lam2_num:.6f}")
print(f"Analytical eigenvalue for k=2:       {analytical_vals[2]:.6f}")
print(f"Relative error:                      {abs(lam2_num - analytical_vals[2]) / analytical_vals[2]:.2e}")
print()

# Check eigenvector matching (accounting for sign ambiguity v and -v)
print("Eigenvector Comparison:")
# For 1st dominant (magnitude |lambda| = 3.732051):
print("1st Dominant Mode (|lambda| = 3.732051):")
print(f"  Numerical v1:               {[round(x, 4) for x in v1_num]}")
print(f"  Analytical v_norm (k=1):     {[round(x, 4) for x in analytical_vecs_norm[1]]}")
print()
# For 2nd dominant (magnitude |lambda| = 3.000000):
print("2nd Dominant Mode (|lambda| = 3.000000):")
print(f"  Numerical v2:               {[round(x, 4) for x in v2_num]}")
print(f"  Analytical v_norm (k=2):     {[round(x, 4) for x in analytical_vecs_norm[2]]}")
print()
print("Note on Sign Discrepancy:")
print("The analytical formula assumes positive coupling between adjacent elements, yielding all positive components for k=1.")
print("However, our matrix has negative off-diagonal entries (-1), which inverts the spatial phase and causes the")
print("highest-energy modes found by Power Iteration to alternate in sign (+, -, +, -, +) across adjacent coordinates.")
print("Additionally, numerical eigenvectors are defined up to a global +/- sign factor.")

