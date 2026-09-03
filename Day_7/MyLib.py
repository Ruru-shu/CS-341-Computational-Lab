# Please refer to comments for explanation on functions and certain variables.
# For Cholesky matrices, instead of storing the whole n*n entries where the upper triangular entries are all 0,
# instead I compress them to a flat list using Cholesky_compress() which can be read via Cholesky_read() and Transpose()
# to get L and L^T.
# This helps store only n(n+1)/2 entries instead of n^2 and saves space. The matrices can be recreated whenever required
# and immediately passed to the required function. They are only stored in memory as such for the required computation,
# and never to disk. The disk always stores a flat list with half the elements.
import os
import subprocess
import math
import sys

class MyComplex:

	def __init__(self, real: float, imag: float):

		self.real = real
		self.imag = imag

	def __repr__(self) -> str:

		sign = "+" if self.imag >= 0 else "-"
		return f"MyComplex({self.real}, {self.imag})"

	def __str__(self) -> str:

		sign = "+" if self.imag >= 0 else "-"
		return f"{self.real} {sign} {abs(self.imag)}j"

	def modulus(self) -> float:

		return round(math.sqrt(self.real**2 + self.imag**2), 3)

	def __add__(self, other: 'MyComplex') -> 'MyComplex':

		return MyComplex(round(self.real + other.real, 3), round(self.imag + other.imag, 3))

	def __sub__(self, other: 'MyComplex') -> 'MyComplex':

		return MyComplex(round(self.real - other.real, 3), round(self.imag - other.imag, 3))

	def __mul__(self, other: 'MyComplex') -> 'MyComplex':

		real_part = round((self.real * other.real) - (self.imag * other.imag), 3)
		imag_part = round((self.real * other.imag) + (self.imag * other.real), 3)
		return MyComplex(real_part, imag_part)

def PRNG(seed=42.0, iter=1):			#Linear Congruential Generator (LCG)

	a = 1103515245.0
	c = 12345.0
	m = 32768.0

	if not hasattr(PRNG, "x0"):
		PRNG.x0 = 42.0
	if seed:
		PRNG.x0 = seed
	for i in range(iter):
		x1 = (a * PRNG.x0 + c) % m
		PRNG.x0 = x1
		yield (x1/m)

def GJE(A, B):     # A is the coefficient matrix; B is the whole number value of the coefficient matrix

	n = len(A)
	for i in range(n):
		A[i].append(B[i])
	for col in range(n):
		max_row = max(range(col, n), key=lambda r: abs(A[r][col]))
		A[col], A[max_row] = A[max_row], A[col]
		pivot = A[col][col]
		for j in range(col, n+1):
			A[col][j] /= pivot
		for row in range(n):
			if row != col:
				factor = A[row][col]
				for j in range(col, n+1):
					A[row][j] -= factor * A[col][j]
	return [A[i][n] for i in range(n)]

def LU(A):		# Doolittle; overwrites A with L and U; L diagonal implicitly 1.0; returns A and perm

	n = len(A)
	perm = list(range(n))
	for i in range(n):
		max_row = max(range(i, n), key=lambda r: abs(A[r][i]))
		A[i], A[max_row] = A[max_row], A[i]
		perm[i], perm[max_row] = perm[max_row], perm[i]
		for j in range(i+1):
			A[j][i] -= sum(A[j][k]*A[k][i] for k in range(j))
		for j in range(i+1, n):
			A[j][i] = (A[j][i] - sum(A[j][k]*A[k][i] for k in range(i))) / A[i][i]
	return A, perm

def forward(A, B):

	if type(A[0]) != list:		# if True, that means Cholesky matrix, and rebuilds the matrix and sets cholesky to True
		A = Cholesky_read(A)
		cholesky = True			# variable that sets if input matrix is Cholesky, otherwise LU
	else:
		cholesky = False
		for row in A:
			if "CHOLESKY" in row:	# Checks for Sentinel string
				cholesky = True
				break
	n = len(A)
	y = [0.0]*n
	for i in range(n):
		y[i] = B[i] - sum(A[i][j]*y[j] for j in range(i))
		if cholesky:
			y[i] /= A[i][i]
	return y

def backward(A, y):

	n = len(A)
	x = [0.0]*n
	for i in range(n-1, -1, -1):
		x[i] = (y[i] - sum(A[i][j]*x[j] for j in range(i+1, n))) / A[i][i]
	return x

def Cholesky(A):	# for HPD matrices; overwrites lower triangle with L, zeros upper

	n = len(A)
	for i in range(n):
		for j in range(n):
			if A[i][j] != A[j][i]:
				sys.exit("Error: matrix is not symmetric")
	for i in range(n):
		A[i][i] = math.sqrt(A[i][i] - sum(A[i][k]**2 for k in range(i)))
		for j in range(i+1, n):
			A[j][i] = (A[j][i] - sum(A[j][k]*A[i][k] for k in range(i))) / A[i][i]
	for i in range(n):
		for j in range(i+1, n):
			A[i][j] = "CHOLESKY"			# Employs a Sentinel String in the upper triangular half, to identify as Cholesky
	return A

def Cholesky_compress(A):		# packs lower triangle (with diag) of a Cholesky matrix into a flat list

	n = len(A)
	flat = []
	for i in range(n):
		for j in range(i+1):
			flat.append(A[i][j])
	return flat

def Cholesky_read(flat, upper=None):		# rebuilds full matrix from a compressed flat list; upper=1 gives L^T instead of L
											# otherwise just using Transpose(Cholesky_read()) gives L^T as well, and as such, 'upper' is redundant. It's only there for readability and uniformity
	n = int((-1 + math.sqrt(1 + 8*len(flat))) / 2)
	A = [[0.0]*n for i in range(n)]
	idx = 0
	for i in range(n):
		for j in range(i+1):
			if upper:
				A[j][i] = flat[idx]
			else:
				A[i][j] = flat[idx]
			idx += 1
	return A

def Pivoting(A, B):	# rearranges rows so no diagonal element of A is 0; B permuted the same way; returns None, None if impossible

	n = len(A)
	A = [row[:] for row in A]
	B = B[:]
	for i in range(n):
		if A[i][i] == 0.0:
			found = False
			for r in range(i+1, n):
				if A[r][i] != 0.0:
					A[i], A[r] = A[r], A[i]
					B[i], B[r] = B[r], B[i]
					found = True
					break
			if not found:
				return None, None
	return A, B

def Jacobi(A, B, tolerance=1e-8):
	A, B = Pivoting(A, B)
	if A is None:
		print("Cannot be pivoted such that no diagonal element is 0")
		return
	n = len(A)
	x = [0.0]*n
	iter = 0
	while True:
		x_new = [0.0]*n
		for i in range(n):
			total = 0.0
			for j in range(n):
				if j != i:
					total += A[i][j] * x[j]
			x_new[i] = (B[i] - total) / A[i][i]

		diff = [x_new[i] - x[i] for i in range(len(x))]		# x_(k+1) - x_k as a vector, diff
		mag = math.sqrt(matmul([diff], Transpose([diff]))[0][0])		# modulus of diff vector

		x = x_new
		iter += 1
		if mag < tolerance:
			break
	print(f"Jacobi converged in {iter} iterations")
	return x

def Gauss_Seidel(A, B, w=1.0, tolerance=1e-8):		# w=1.0 gives plain Gauss-Seidel; any other w gives SOR
	A, B = Pivoting(A, B)
	if A is None:
		print("Cannot be pivoted such that no diagonal element is 0")
		return
	n = len(A)
	x = [0.0]*n
	iter = 0
	while True:
		temp = []							# holds x_(k+1)[i] - x_k[i] for this iter only
		for i in range(n):
			total = 0.0
			for j in range(n):
				if j != i:
					total += A[i][j] * x[j]
			new_val = (B[i] - total) / A[i][i]		# plain Gauss-Seidel value
			old_val = x[i]
			final_val = (1-w)*old_val + w*new_val	# SOR blend of old and new; final_val = new_val when w=1.0
			temp.append(final_val - old_val)
			x[i] = final_val					# overwrite immediately, same array

		mag = math.sqrt(matmul([temp], Transpose([temp]))[0][0])
		del temp
		iter += 1
		if mag < tolerance:
			break
	print(f"Gauss-Seidel (w={w}) converged in {iter} iterations")
	return x

def Bracketing(f, a, b, beta=0.5):		# finds proper bracketting and is called at the start of Regula_Falsi() and Bisection(); and returns proper [a, b]
	while f(a)*f(b) > 0:				# shifts a or b outward until f(a)*f(b) < 0
		if abs(f(a)) < abs(f(b)):
			a = a - beta*(b-a)
		else:
			b = b + beta*(b-a)
	return a, b

def Bisection(f, a, b, beta=0.5, tolerance=1e-6):
	if f(a)*f(b) > 0:
		a, b = Bracketing(f, a, b, beta)		# To get proper [a.b]
	iters = 0
	while (b-a) > tolerance:
		c = (a+b)/2
		if f(c)*f(a) < 0:
			b = c
		else:
			a = c
		iters += 1
	print(f"Bisection converged in {iters} iterations")
	return a if abs(f(a)) < abs(f(b)) else b

def Regula_Falsi(f, a, b, beta=0.5, tolerance=1e-6):
	if f(a)*f(b) > 0:
		a, b = Bracketing(f, a, b, beta)	# To get proper [a.b]
	iters = 0
	c_prev = a
	while True:
		c = b - (b-a)*f(b)/(f(b)-f(a))
		if f(a)*f(c) < 0:
			b = c
		else:
			a = c
		iters += 1
		if abs(c - c_prev) < tolerance:
			break
		c_prev = c
	print(f"Regula Falsi converged in {iters} iterations")
	return a if abs(f(a)) < abs(f(b)) else b

def Transpose(A):			#Transposes a given matrix

	rows = len(A)
	cols = len(A[0])
	return [[A[i][j] for i in range(rows)] for j in range(cols)]

def Tri_Mats(A):		# splits A into Lower and Upper triangular (both with diag) as flat lists, Cholesky_read() compatible

	n = len(A)
	lower = []
	upper = []
	for i in range(n):
		for j in range(i+1):
			lower.append(A[i][j])
			upper.append(A[j][i])
	return lower, upper

def invert(A):			# Inverts a given matrix via GJE

	n = len(A)
	aug = [A[i][:] + [1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
	for col in range(n):
		max_row = max(range(col, n), key=lambda r: abs(aug[r][col]))
		aug[col], aug[max_row] = aug[max_row], aug[col]
		pivot = aug[col][col]
		for j in range(col, 2*n):
			aug[col][j] /= pivot
		for row in range(n):
			if row != col:
				factor = aug[row][col]
				for j in range(col, 2*n):
					aug[row][j] -= factor * aug[col][j]
	return [aug[i][n:] for i in range(n)]

def Determinant(A):

	if type(A[0]) != list:		# compressed Cholesky flat list, diagonals sit at index 0, 2, 5, 9...
		n = int((-1 + math.sqrt(1 + 8*len(A))) / 2)
		idx = 0
		det = 1.0
		for i in range(n):
			idx += i
			det *= A[idx]
			idx += 1
		return round(det**2, 3)
	cholesky = False
	for row in A:
		if "CHOLESKY" in row:
			cholesky = True
			break
	n = len(A)
	det = 1.0
	for i in range(n):
		det *= A[i][i]
	if cholesky:
		det = det**2
	return round(det, 3)

def matmul(A, B):			# Multiplies 2 given matrices

	rows_A, cols_A = len(A), len(A[0])
	rows_B, cols_B = len(B), len(B[0])
	if cols_A != rows_B:
		sys.exit(f"Error: cannot multiply {rows_A}x{cols_A} matrix with {rows_B}x{cols_B} matrix")
	C = [[0.0]*cols_B for i in range(rows_A)]
	for i in range(rows_A):
		for j in range(cols_B):
			for k in range(cols_A):
				C[i][j] += A[i][k]*B[k][j]
	return C

if __name__ == "__main__":			# Outputs all .py files in the same subfolder to .txt

	folder = os.path.dirname(os.path.abspath(__file__))
	for fname in os.listdir(folder):
		if fname.endswith(".py") and fname != os.path.basename(__file__):
			outname = fname.replace(".py", ".txt")
			outpath = os.path.join(folder, outname)
			with open(outpath, "w") as f:
				subprocess.run(["python", os.path.join(folder, fname)], stdout=f, stderr=f)
			print(f"ran {fname} -> {outname}")
