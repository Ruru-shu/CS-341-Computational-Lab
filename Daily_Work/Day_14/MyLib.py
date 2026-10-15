# Please refer to comments for explanation on functions and certain variables.
# For Cholesky matrices, instead of storing the whole n*n entries where the upper triangular entries are all 0,
# instead I compress them to a flat list using Cholesky_compress() which can be read via Cholesky_read() and Transpose()
# to get L and L^T.
# This helps store only n(n+1)/2 entries instead of n^2 and saves space. The matrices can be recreated whenever required
# and immediately passed to the required function. They are only stored in memory as such for the required computation,
# and never to disk. The disk always stores a flat list with half the elements.
# The Save_Data() and Snapshot() and Report_Files() are helpers to __main__ branch of this file. The __main__ branch when run,
# executes every other .py file in the folder, and saves the text output as .txt files, any 2D plots with matplotlib as .png and any 3D plots as .dat (GNUplot readable)
import os
import subprocess
import math
import sys
import matplotlib
import matplotlib.pyplot as plt

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

def PRNG(seed=None, iter=1):                       # Linear Congruential Generator (LCG)

    # The function attribute PRNG.x0 stores the last generated value.
    # If no seed is supplied, the previous PRNG.x0 is automatically used as the seed, allowing successive calls to continue the same sequence.

    if not hasattr(PRNG, "x0"):
        PRNG.x0 = 42                                # Initial default seed

    if seed is not None:
        PRNG.x0 = seed                              # Use explicitly supplied seed

    a = 1103515245
    c = 12345
    m = 32768
    numbers = []

    for i in range(iter):
        x1 = (a * PRNG.x0 + c) % m
        PRNG.x0 = x1                               # Store last generated value
        numbers.append(x1 / m)

    if iter == 1:
        return numbers[0]                           # Return a float
    else:
        return numbers                              # Return a list

def GJE(A, B):     # A is the coefficient matrix; B is the whole number value of the coefficient matrix
	# B is appended to A as an extra column, forming the Augmented Matrix [A|B], which is then reduced in place.
	# As such, A itself gets overwritten (it ends up as the Identity matrix, with the solution as the extra last column).

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

def Jacobi(A, B, tolerance=1e-8, max_iter=10000):		# max_iter is needed: if A is not diagonally dominant the iteration can diverge, and without a limit it would never stop
	A, B = Pivoting(A, B)
	if A is None:
		print("Cannot be pivoted such that no diagonal element is 0")
		return
	n = len(A)
	x = [0.0]*n
	iter = 0
	converged = False
	while iter < max_iter:
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
			converged = True
			break
	if converged:
		print(f"Jacobi converged in {iter} iterations")
	else:
		print(f"Jacobi did not converge in {max_iter} iterations")
	return x

def Gauss_Seidel(A, B, w=1.0, tolerance=1e-8, max_iter=10000):		# w=1.0 gives plain Gauss-Seidel; any other w gives SOR
	# max_iter is needed: if A is not diagonally dominant the iteration can diverge, and without a limit it would never stop
	A, B = Pivoting(A, B)
	if A is None:
		print("Cannot be pivoted such that no diagonal element is 0")
		return
	n = len(A)
	x = [0.0]*n
	iter = 0
	converged = False
	while iter < max_iter:
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
			converged = True
			break
	if converged:
		print(f"Gauss-Seidel (w={w}) converged in {iter} iterations")
	else:
		print(f"Gauss-Seidel (w={w}) did not converge in {max_iter} iterations")
	return x

def Bracketing(f, a, b, beta=0.5, max_iter=100):		# finds proper bracketting and is called at the start of Regula_Falsi() and Bisection(); and returns proper [a, b]
	# max_iter is needed: if f never changes sign (e.g. x^2 + 1) the interval would be widened forever. 100 widenings already stretch it by a factor of ~1e17
	# returns None, None if no bracket was found within max_iter
	iters = 0
	while f(a)*f(b) > 0 and iters < max_iter:				# shifts a or b outward until f(a)*f(b) < 0
		if abs(f(a)) < abs(f(b)):
			a = a - beta*(b-a)
		else:
			b = b + beta*(b-a)
		iters += 1
	if f(a)*f(b) > 0:
		return None, None
	return a, b

def Bisection(f, a, b, beta=0.5, tolerance=1e-8, max_iter=10000):		# max_iter stops the loop if the tolerance can never be reached (e.g. a tolerance below the float precision)
	if f(a)*f(b) > 0:
		a, b = Bracketing(f, a, b, beta)		# To get proper [a.b]
	if a is None:
		print("Could not find a bracketing interval")
		return
	iters = 0
	while (b-a) > tolerance and iters < max_iter:
		c = (a+b)/2
		if f(c)*f(a) < 0:
			b = c
		else:
			a = c
		iters += 1
	if (b-a) > tolerance:
		print(f"Bisection did not converge in {max_iter} iterations")
	else:
		print(f"Bisection converged in {iters} iterations")
	return a if abs(f(a)) < abs(f(b)) else b

def Regula_Falsi(f, a, b, beta=0.5, tolerance=1e-8, max_iter=10000):		# max_iter stops the loop if the tolerance can never be reached (e.g. a tolerance below the float precision)
	if f(a)*f(b) > 0:
		a, b = Bracketing(f, a, b, beta)	# To get proper [a.b]
	if a is None:
		print("Could not find a bracketing interval")
		return
	iters = 0
	c_prev = a
	converged = False
	while iters < max_iter:
		c = b - (b-a)*f(b)/(f(b)-f(a))
		if f(a)*f(c) < 0:
			b = c
		else:
			a = c
		iters += 1
		if abs(c - c_prev) < tolerance:
			converged = True
			break
		c_prev = c
	if converged:
		print(f"Regula Falsi converged in {iters} iterations")
	else:
		print(f"Regula Falsi did not converge in {max_iter} iterations")
	return a if abs(f(a)) < abs(f(b)) else b

def Fixed_Point(g, x0, tolerance=1e-8, max_iter=10000):		# g(x) is the rearranged form of f(x) = 0 as x = g(x); x0 is initial guess
	iters = 0
	x = x0
	while iters < max_iter:
		x_new = g(x)
		iters += 1
		if abs(x_new - x) < tolerance:
			break
		x = x_new
	if iters == max_iter:
		print(f"Fixed Point did not converge in {max_iter} iterations")
	else:
		print(f"Fixed Point converged in {iters} iterations")
	return x_new

def Newton_Raphson(F, J, x0, tolerance=1e-8, max_iter=10000, verbose=True):			# it takes for both single variable and multivariable.
	# F is list of functions for multivariable; or single function for 1D
	# J is 2D Jacobian matrix of functions; or single derivative function for 1D
	# x0 is the initial guess.
	# x0 is used to determine if it's 1 var or multivar. x0 is a float if single var, and a list of floats if multivar.
	# verbose: 'True' prints the "converged / did not converge" line.
	#          'False' suppresses that print. It is meant for when Newton_Raphson() is called many times from inside another
	#          function (Euler_Bkwd()), where hundreds of identical lines would flood the output.
	iters = 0
	if type(x0) is not list:				# 1D case
		x = x0
		while iters < max_iter:
			fx = F(x)
			jx = J(x)
			x_new = x - fx/jx
			iters += 1
			if abs(x_new - x) < tolerance:
				break
			x = x_new
		if verbose:
			if iters == max_iter:
				print(f"Newton-Raphson did not converge in {max_iter} iterations")
			else:
				print(f"Newton-Raphson converged in {iters} iterations")
		return x_new
	else:									# multivariable case
		n = len(x0)
		x = x0[:]
		while iters < max_iter:
			fx = [F[i](x) for i in range(n)]				# evaluates F at current x
			Jx = [[J[i][j](x) for j in range(n)] for i in range(n)]	# evaluates Jacobian at current x
			Jinv = invert(Jx)									# inverts J(x)
			dx = [sum(Jinv[i][j]*fx[j] for j in range(n)) for i in range(n)]	# x_new = x - J^(-1)*x*f(x)
			x_new = [x[i] - dx[i] for i in range(n)]
			iters += 1
			mag = math.sqrt(sum((x_new[i]-x[i])**2 for i in range(n)))		# norm of x, y over old x, y difference
			if mag < tolerance:
				break
			x = x_new
		if verbose:
			if iters == max_iter:
				print(f"Newton-Raphson did not converge in {max_iter} iterations")
			else:
				print(f"Newton-Raphson converged in {iters} iterations")
		return x_new

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

def Determinant(A, perm=None, raw=False):
	# Cholesky input (flat list, or matrix carrying the "CHOLESKY" sentinel): det = (product of the diagonal of L)^2
	# Otherwise A is taken to be an already decomposed (LU) matrix: det = product of the diagonal of U, times the sign of perm
	# perm is the list returned by LU(); every row swap LU() made flips the sign. perm=None means no swaps are accounted for
	# raw=True: A is an ordinary matrix, and the determinant is found on a copy of it (A is not modified)
	# The result is rounded to 12 significant digits, which hides float noise (23.999999999999996 -> 24.0) without wiping out small determinants

	if type(A[0]) != list:		# compressed Cholesky flat list, diagonals sit at index 0, 2, 5, 9...
		n = int((-1 + math.sqrt(1 + 8*len(A))) / 2)
		idx = 0
		det = 1.0
		for i in range(n):
			idx += i
			det *= A[idx]
			idx += 1
		return float(f"{det**2:.12g}")
	cholesky = False
	for row in A:
		if "CHOLESKY" in row:
			cholesky = True
			break
	n = len(A)
	if raw and not cholesky:
		M = [row[:] for row in A]
		det = 1.0
		for col in range(n):
			max_row = col
			for r in range(col+1, n):
				if abs(M[r][col]) > abs(M[max_row][col]):
					max_row = r
			if M[max_row][col] == 0:
				return 0.0
			if max_row != col:
				M[col], M[max_row] = M[max_row], M[col]
				det = -det
			det *= M[col][col]
			for r in range(col+1, n):
				factor = M[r][col] / M[col][col]
				for j in range(col, n):
					M[r][j] -= factor * M[col][j]
		return float(f"{det:.12g}")
	det = 1.0
	for i in range(n):
		det *= A[i][i]
	if cholesky:
		det = det**2
	if perm is not None and not cholesky:
		for i in range(n):		# sign of perm: every pair that is out of order counts as one swap
			for j in range(i+1, n):
				if perm[i] > perm[j]:
					det = -det
	return float(f"{det:.12g}")

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

def deflate(coeff, root):		# synthetic division; removes found root from polynomial; returns reduced coefficients
	new_coeff = [coeff[0]]
	for i in range(1, len(coeff)):
		new_coeff.append(coeff[i] + root*new_coeff[-1])
	return new_coeff[:-1]		# last entry is remainder (should be ~0 if root is correct)

def derivative(coeff):			# returns coefficients of derivative of a polynomial
	n = len(coeff) - 1
	dcoeff = []
	for i in range(n):
		dcoeff.append(coeff[i] * (n-i))
	return dcoeff

def poly_eval(coeff, x):		# evaluates polynomial at x; coefficients in descending order of power
	result = 0.0
	for i in range(len(coeff)):
		result += coeff[i] * (x ** (len(coeff)-1-i))
	return result

def Laguerre(coeff, x0=0.0, tolerance=1e-8, max_iter=10000):
	roots = []
	n_orig = len(coeff) - 1
	while len(coeff) > 1:
		n = len(coeff) - 1
		x = x0
		iters = 0
		while iters < max_iter:
			fx = poly_eval(coeff, x)
			if abs(fx) < tolerance:
				break
			d1 = derivative(coeff)
			d2 = derivative(d1)
			g = poly_eval(d1, x) / fx
			h = g**2 - poly_eval(d2, x) / fx
			disc = math.sqrt(abs((n-1) * (n*h - g**2)))		#discriminant
			denom_plus  = g + disc
			denom_minus = g - disc
			if abs(denom_plus) >= abs(denom_minus):
				a = n / denom_plus
			else:
				a = n / denom_minus
			x_new = x - a
			iters += 1
			if abs(x_new - x) < tolerance:
				x = x_new
				break
			x = x_new

		if iters == max_iter:
			print(f"Laguerre did not converge for this root")
		roots.append(round(x, 6))
		coeff = deflate(coeff, x)
	return roots

def poly_str(coeff):			# Shows a polynomial neatly
    n = len(coeff) - 1
    terms = []
    for i in range(len(coeff)):
        terms.append(f"{coeff[i]}x^{n-i}")
    return " + ".join(terms)

def Midpoint_Int(f, a, b, N=10000):
	N = int(N)
	h = float((b-a)/N)
	I = 0
	for i in range(N):
		I += f(a + (i + 0.5)*h)*h
	return I

def Trapezoidal_Int(f, a, b, N=10000):
	N = int(N)
	h = float((b-a)/N)
	I = 0.5*(f(a) + f(b))
	for i in range(1, N):
		I += f(a + i*h)
	return I*h

def Simpson_Int(f, a, b, N=10000):
	N = int(N)
	N = N if N % 2 == 0 else N + 1		# Simpson requires even N
	h = (b-a)/N
	s = f(a)+f(b)
	for i in range(1, N):
		x = a + i*h
		if i % 2 == 0:
			s += 2*f(x)
		else:
			s += 4*f(x)
	return (h/3)*s

def MonteCarlo_Int(f, a, b, N=10000):
	N = int(N)
	s = 0
	for i in range(N):
		x = PRNG()*(b-a) + a
		s+= f(x)
	return (b-a)*s/N

def Gaussian_Int(f, a, b, n=5):		# Gauss-Legendre for finite [a,b]; Gauss-Laguerre when b=float('inf')
	# n is the number of nodes, 1 to 10
	legendre = {		# Legendre nodes and weights on [-1, 1] for n = 1 to 10
		1: ([0.0],
			[2.0]),
		2: ([-0.5773502691896257, 0.5773502691896257],
			[1.0, 1.0]),
		3: ([-0.7745966692414834, 0.0, 0.7745966692414834],
			[0.5555555555555556, 0.8888888888888888, 0.5555555555555556]),
		4: ([-0.8611363115940526, -0.33998104358485626, 0.33998104358485626, 0.8611363115940526],
			[0.34785484513745385, 0.6521451548625461, 0.6521451548625461, 0.34785484513745385]),
		5: ([-0.906179845938664, -0.5384693101056831, 0.0, 0.5384693101056831, 0.906179845938664],
			[0.23692688505618908, 0.47862867049936647, 0.5688888888888889, 0.47862867049936647, 0.23692688505618908]),
		6: ([-0.932469514203152, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, 0.6612093864662645,
				0.932469514203152],
			[0.17132449237917036, 0.3607615730481386, 0.46791393457269104, 0.46791393457269104, 0.3607615730481386,
				0.17132449237917036]),
		7: ([-0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0.0, 0.4058451513773972,
				0.7415311855993945, 0.9491079123427585],
			[0.1294849661688697, 0.27970539148927664, 0.3818300505051189, 0.4179591836734694, 0.3818300505051189,
				0.27970539148927664, 0.1294849661688697]),
		8: ([-0.9602898564975363, -0.7966664774136267, -0.525532409916329, -0.1834346424956498, 0.1834346424956498,
				0.525532409916329, 0.7966664774136267, 0.9602898564975363],
			[0.10122853629037626, 0.22238103445337448, 0.31370664587788727, 0.362683783378362, 0.362683783378362,
				0.31370664587788727, 0.22238103445337448, 0.10122853629037626]),
		9: ([-0.9681602395076261, -0.8360311073266358, -0.6133714327005904, -0.3242534234038089, 0.0,
				0.3242534234038089, 0.6133714327005904, 0.8360311073266358, 0.9681602395076261],
			[0.08127438836157441, 0.1806481606948574, 0.26061069640293544, 0.31234707704000286, 0.3302393550012598,
				0.31234707704000286, 0.26061069640293544, 0.1806481606948574, 0.08127438836157441]),
		10: ([-0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472, -0.14887433898163122,
				0.14887433898163122, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845, 0.9739065285171717],
			[0.06667134430868814, 0.1494513491505806, 0.21908636251598204, 0.26926671930999635, 0.29552422471475287,
				0.29552422471475287, 0.26926671930999635, 0.21908636251598204, 0.1494513491505806, 0.06667134430868814])
	}

	laguerre = {		# Laguerre nodes and weights for n = 1 to 10
		1: ([1.0],
			[1.0]),
		2: ([0.585786437626905, 3.414213562373095],
			[0.8535533905932737, 0.14644660940672624]),
		3: ([0.4157745567834791, 2.294280360279042, 6.2899450829374794],
			[0.711093009929173, 0.27851773356924087, 0.010389256501586135]),
		4: ([0.3225476896193923, 1.7457611011583465, 4.536620296921128, 9.395070912301133],
			[0.6031541043416336, 0.35741869243779967, 0.038887908515005384, 0.0005392947055613274]),
		5: ([0.2635603197181409, 1.4134030591065168, 3.596425771040722, 7.085810005858837, 12.640800844275782],
			[0.5217556105828086, 0.3986668110831759, 0.0759424496817076, 0.0036117586799220484, 2.3369972385776228e-05]),
		6: ([0.2228466041792607, 1.188932101672623, 2.992736326059314, 5.77514356910451, 9.83746741838259,
				15.982873980601703],
			[0.4589646739499636, 0.41700083077212097, 0.11337338207404497, 0.010399197453149074, 0.00026101720281493206,
				8.985479064296212e-07]),
		7: ([0.19304367656036242, 1.026664895339192, 2.567876744950746, 4.900353084526484, 8.18215344456286,
				12.734180291797815, 19.39572786226254],
			[0.4093189517012739, 0.42183127786171976, 0.1471263486575053, 0.02063351446871694, 0.0010740101432807456,
				1.5865464348564203e-05, 3.17031547899558e-08]),
		8: ([0.170279632305101, 0.9037017767993799, 2.2510866298661307, 4.266700170287659, 7.0459054023934655,
				10.758516010180996, 15.740678641278004, 22.863131736889265],
			[0.3691885893416375, 0.41878678081434295, 0.1757949866371718, 0.03334349226121565, 0.0027945362352256725,
				9.076508773358213e-05, 8.485746716272531e-07, 1.0480011748715104e-09]),
		9: ([0.15232222773180826, 0.8072200227422558, 2.0051351556193473, 3.783473973331233, 6.204956777876613,
				9.372985251687576, 13.466236911092093, 18.833597788991696, 26.374071890927375],
			[0.33612642179796254, 0.4112139804239844, 0.1992875253708856, 0.0474605627656516, 0.005599626610794583,
				0.00030524976709321057, 6.592123026075352e-06, 4.110769330349548e-08, 3.290874030350708e-11]),
		10: ([0.13779347054049243, 0.7294545495031705, 1.808342901740316, 3.4014336978548996, 5.552496140063804,
				8.330152746764497, 11.843785837900066, 16.279257831378104, 21.99658581198076, 29.92069701227389],
			[0.30844111576502015, 0.40111992915527356, 0.2180682876118094, 0.062087456098677746, 0.0095015169751811,
				0.0007530083885875388, 2.8259233495995656e-05, 4.2493139849626863e-07, 1.8395648239796308e-09, 9.911827219609008e-13])
	}

	if n < 1 or n > 10:
		sys.exit("Error: n must be between 1 and 10")
	if b == float('inf'):				# Gauss-Laguerre branch; f is passed as e^-x * g(x), so we extract g(x) = f(x)*e^x
		nodes, weights = laguerre[n]
		return sum(weights[i] * f(nodes[i] + a) * math.exp(nodes[i]) for i in range(n))
									# the integral from a to inf is shifted to 0 to inf, hence f(node + a); e^-x weight is built into the quadrature; we undo e^-x from f to get bare g(x)
	else:								# Gauss-Legendre branch; change of variable from [a,b] to [-1,1]
		nodes, weights = legendre[n]
		mid   = (b + a) / 2.0			# (b+a)/2 shift
		half  = (b - a) / 2.0			# (b-a)/2 scale
		return half * sum(weights[i] * f(half * nodes[i] + mid) for i in range(n))

def Euler_Fwd(f, x0, y0, x_end, h, tolerance=None):		# Forward (explicit) Euler for dy/dx = f(x, y); returns lists xs, ys
	# tolerance=None (default): plain Forward Euler. If given, each step is redone with 2, 4, 8... sub-steps until the result stops changing by tolerance
	xs = [x0]
	ys = [y0]
	x = x0
	y = y0
	while x < x_end - 1e-12:
		step = min(h, x_end - x)		# last step is clipped to land exactly on x_end
		y_new = y + step*f(x, y)
		if tolerance is not None:
			m = 1
			change = tolerance + 1
			while change >= tolerance and m < 1024:
				m = 2*m
				x_try = x
				y_try = y
				for k in range(m):
					y_try += (step/m)*f(x_try, y_try)
					x_try += step/m
				change = abs(y_try - y_new)
				y_new = y_try
		x += step
		y = y_new
		xs.append(x)
		ys.append(y)
	return xs, ys

def Euler_Bkwd(f, x0, y0, x_end, h, tolerance=1e-8):		# Backward (implicit) Euler for dy/dx = f(x, y); returns lists xs, ys
	# each step solves y_(n+1) = y_n + h*f(x_(n+1), y_(n+1)) using Newton_Raphson(); tolerance is passed on to it
	xs = [x0]
	ys = [y0]
	x = x0
	y = y0
	while x < x_end - 1e-12:
		step = min(h, x_end - x)
		x1 = x + step
		def G(z):
			return z - y - step*f(x1, z)
		def dG(z):
			return 1 - step*(f(x1, z + 1e-6) - f(x1, z - 1e-6))/2e-6
		z = Newton_Raphson(G, dG, y + step*f(x, y), tolerance, 100, False)
		if abs(G(z)) > 1e-6:
			break
		x = x1
		y = z
		xs.append(x)
		ys.append(y)
	return xs, ys

def Predictor_Corrector(f, x0, y0, x_end, h, tolerance=None):		# Euler predictor + trapezoidal corrector; returns lists xs, ys
	# tolerance=None (default): corrector applied once per step. If given, corrector is repeated until it changes by less than tolerance
	xs = [x0]
	ys = [y0]
	x = x0
	y = y0
	while x < x_end - 1e-12:
		step = min(h, x_end - x)
		x1 = x + step
		slope = f(x, y)
		y_pred = y + step*slope							# predictor (Forward Euler)
		y_new = y + 0.5*step*(slope + f(x1, y_pred))		# corrector (trapezoidal)
		if tolerance is not None:
			for i in range(100):
				y_next = y + 0.5*step*(slope + f(x1, y_new))
				change = abs(y_next - y_new)
				y_new = y_next
				if change < tolerance:
					break
		x = x1
		y = y_new
		xs.append(x)
		ys.append(y)
	return xs, ys

def RK4(f, x0, Y0, x_end, h):		# Runge-Kutta 4 for dY/dx = f(x, Y); returns lists xs, ys
	# Y0 is a float for a 1st order ODE (f returns a float), or a list for a 2nd order ODE / system (f returns a list of slopes, one per variable)
	# For a list Y0, ys is a list of lists: ys[n] holds all the variables at xs[n]
	xs = [x0]
	ys = [Y0]
	x = x0
	Y = Y0
	system = type(Y0) == list
	while x < x_end - 1e-12:
		step = min(h, x_end - x)		# last step is clipped to land exactly on x_end
		if system:
			n = len(Y)
			k1 = f(x, Y)
			k2 = f(x + step/2, [Y[i] + step/2*k1[i] for i in range(n)])
			k3 = f(x + step/2, [Y[i] + step/2*k2[i] for i in range(n)])
			k4 = f(x + step, [Y[i] + step*k3[i] for i in range(n)])
			Y = [Y[i] + step/6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(n)]
		else:
			k1 = f(x, Y)
			k2 = f(x + step/2, Y + step/2*k1)
			k3 = f(x + step/2, Y + step/2*k2)
			k4 = f(x + step, Y + step*k3)
			Y = Y + step/6*(k1 + 2*k2 + 2*k3 + k4)
		x += step
		xs.append(x)
		ys.append(Y)
	return xs, ys

def Shooting(f, a, b, ya, yb, s1, s2, h, tolerance=1e-8, max_iter=100):		# Shooting method for a 2nd order BVP y'' = g(x, y, y'), y(a) = ya, y(b) = yb
	# f(x, Y) is the same system form that RK4() takes: Y = [y, y'] and f returns [y', y'']
	# s1, s2 are the two initial guesses for the slope y'(a); h is the RK4 step
	# every shot is solved with RK4(). s1, s2 must bracket yb: one guess undershoots (y(b) < yb), the other overshoots (y(b) > yb); otherwise returns None, None
	# the next slope comes from Lagrange's linear interpolation between the low guess zl and the high guess zh:
	#     zeta = zl + (zh - zl)/(yh - yl) * (yb - yl)
	# the new shot then replaces zl if it undershoots, or zh if it overshoots, so the bracket always closes in
	# returns xs, shots: shots is a list with one entry per shot (each entry is the ys list of RK4()), so every trajectory can be plotted
	# the last shot is the converged solution; its slope is shots[-1][0][1]
	shots = []
	xs, Y1 = RK4(f, a, [ya, s1], b, h)
	shots.append(Y1)
	xs, Y2 = RK4(f, a, [ya, s2], b, h)
	shots.append(Y2)
	end1 = Y1[-1][0]
	end2 = Y2[-1][0]
	if (end1 - yb)*(end2 - yb) > 0:
		print("Shooting: the two guesses do not bracket the boundary value")
		return None, None
	if end1 < end2:
		zl, yl, zh, yh = s1, end1, s2, end2
	else:
		zl, yl, zh, yh = s2, end2, s1, end1
	iters = 0
	r = tolerance + 1
	while abs(r) > tolerance and iters < max_iter:
		zeta = zl + (zh - zl)/(yh - yl)*(yb - yl)		# Lagrange linear interpolation for the next slope
		xs, Y3 = RK4(f, a, [ya, zeta], b, h)
		shots.append(Y3)
		y_end = Y3[-1][0]
		r = y_end - yb
		if y_end < yb:
			zl = zeta
			yl = y_end
		else:
			zh = zeta
			yh = y_end
		iters += 1
	if abs(r) > tolerance:
		print(f"Shooting did not converge in {max_iter} iterations")
	else:
		print(f"Shooting converged in {len(shots)} shots")
	return xs, shots

def Finite_Elim(f, a, b, ya, yb, N):		# Finite element method as in the lecture, for y'' = f(x), y(a) = ya, y(b) = yb
	# N is the number of equal sub-intervals, so there are N+1 nodes; returns the nodes xs and the nodal values ys
	# y''(x_n) ~ (y_(n+1) - 2 y_n + y_(n-1))/h^2 gives, for the N-1 interior nodes, the tridiagonal symmetric positive definite system
	#     -y_(n-1) + 2 y_n - y_(n+1) = -h^2 f(x_n)
	# ya and yb are moved to the right side of the first and last equation
	# the system is solved by Cholesky(); the factor is packed with Cholesky_compress(), forward() rebuilds L from that flat list,
	# and backward() is given L^T from Cholesky_read(upper=1)
	h = (b-a)/N
	xs = [a + i*h for i in range(N+1)]
	n = N - 1
	A = [[0.0]*n for i in range(n)]
	B = [0.0]*n
	for i in range(n):
		A[i][i] = 2.0
		if i > 0:
			A[i][i-1] = -1.0
		if i < n-1:
			A[i][i+1] = -1.0
		B[i] = -h*h*f(xs[i+1])
	B[0] += ya
	B[n-1] += yb
	Cholesky(A)
	flat = Cholesky_compress(A)
	y = forward(flat, B)
	inner = backward(Cholesky_read(flat, upper=1), y)
	ys = [ya] + inner + [yb]
	return xs, ys

def Save_Data(filename, *columns):		# writes the given equal-length lists as whitespace separated columns to <filename> in the same folder as MyLib.py
	# readable by gnuplot, e.g. splot 'Q4.dat' using 1:2:3 with lines
	folder = os.path.dirname(os.path.abspath(__file__))
	with open(os.path.join(folder, filename), "w") as f:
		for i in range(len(columns[0])):
			row = [str(col[i]) for col in columns]
			f.write(" ".join(row) + "\n")

def Snapshot(folder):		# returns {file name: last modified time} for every file in the folder; two snapshots compared show which files were created or changed
	snap = {}
	for name in os.listdir(folder):
		path = os.path.join(folder, name)
		if os.path.isfile(path):
			snap[name] = os.path.getmtime(path)
	return snap

def Report_Files(before, after, skip=None):		# prints every file that is new or changed between two Snapshot()s; skip is a file name to leave out (the .txt, which is already announced)
	for name in sorted(after):
		if name == skip:
			continue
		if name not in before:
			print(f"   created {name}")
		elif after[name] != before[name]:
			print(f"   updated {name}")

if __name__ == "__main__":			# Outputs all .py files in the same subfolder to .txt, and the matplotlib figures they make to .png (3D plots also to .dat)

	folder = os.path.dirname(os.path.abspath(__file__))
	if len(sys.argv) > 1:			# this branch is entered by the loop below, with one script as argument
		import runpy				# runs that script with the Agg backend (no plot windows), then saves all its figures as <script>_1.png, <script>_2.png...
		import warnings
		import matplotlib
		matplotlib.use("Agg")
		import matplotlib.pyplot as plt
		warnings.filterwarnings("ignore", message=".*non-interactive")		# plt.show() does nothing under Agg
		quiet = len(sys.argv) > 2		# the batch loop below passes a 2nd argument, so that the file lines are not captured into the .txt; it reports the files itself
		script_folder = os.path.dirname(os.path.abspath(sys.argv[1]))		# the figures are saved next to the script
		before = Snapshot(script_folder)
		runpy.run_path(sys.argv[1], run_name="__main__")
		name = os.path.splitext(sys.argv[1])[0]
		k = 1
		for num in plt.get_fignums():
			fig = plt.figure(num)
			fig.savefig(f"{name}_{k}.png")
			for ax in fig.axes:				# 3D axes: also saves the plotted line(s) as <script>_<k>.dat, for gnuplot (splot '<file>.dat' using 1:2:3 with lines)
				if hasattr(ax, "get_zlim"):
					for j in range(len(ax.lines)):
						x, y, z = ax.lines[j].get_data_3d()
						suffix = "" if len(ax.lines) == 1 else f"_{j+1}"
						Save_Data(f"{name}_{k}{suffix}.dat", x, y, z)
			k += 1
		if not quiet:
			Report_Files(before, Snapshot(script_folder))
	else:
		for fname in os.listdir(folder):
			if fname.endswith(".py") and fname != os.path.basename(__file__):
				outname = fname.replace(".py", ".txt")
				outpath = os.path.join(folder, outname)
				before = Snapshot(folder)
				result = subprocess.run(["python", os.path.abspath(__file__), os.path.join(folder, fname), "batch"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
				if result.stdout != "":			# .txt only for scripts that actually print something
					with open(outpath, "w") as f:
						f.write(result.stdout)
					print(f"ran {fname} -> {outname}")
				else:
					print(f"ran {fname}")
				Report_Files(before, Snapshot(folder), outname)
