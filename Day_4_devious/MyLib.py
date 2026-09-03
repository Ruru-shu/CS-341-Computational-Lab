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

def PRNG(seed=42.0, iter=1):
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

def forward(A, B, unit_diag=True):		# unit_diag=True for Doolittle L; False for Cholesky L;
	n = len(A)
	y = [0.0]*n
	for i in range(n):
		y[i] = B[i] - sum(A[i][j]*y[j] for j in range(i))
		if not unit_diag:
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
			A[i][j] = 0.0
	return A

def invert(A):
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

def matmul(A, B):
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

if __name__ == "__main__":
	folder = os.path.dirname(os.path.abspath(__file__))
	for fname in os.listdir(folder):
		if fname.endswith(".py") and fname != os.path.basename(__file__):
			outname = fname.replace(".py", ".txt")
			outpath = os.path.join(folder, outname)
			with open(outpath, "w") as f:
				subprocess.run(["python", os.path.join(folder, fname)], stdout=f, stderr=f)
			print(f"ran {fname} -> {outname}")
