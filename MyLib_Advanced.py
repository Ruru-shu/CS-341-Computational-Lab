# MyLib_Advanced.py
# Extends MyLib.py with derivative-free variants of Newton_Raphson() and Fixed_Point().
# Instead of requiring an analytic Jacobian/derivative (Newton_Raphson) or a rearranged
# g(x) (Fixed_Point), both now take only F (or f) directly, and use Deriv() to estimate
# derivatives numerically via the Secant Method (a two-point finite-difference slope)
# at a small step h.
#
# Note on Deriv(): it only ever returns a single number (the derivative value at a point),
# never a callable. So it can't be recursed into directly to get higher-order derivatives.
# You COULD wrap it in a lambda (e.g. Deriv(lambda t: Deriv(f, t, h), x, h)) to approximate
# f''(x), but this compounds finite-difference error badly (you're differencing differences),
# so it isn't done here.
#
# Fixed_Point() below is deliberately NOT full Newton-Raphson in disguise: for the
# multivariable case it only ever uses each F[i]'s own self-derivative (dF_i/dx_i), the
# same way Jacobi() in MyLib.py only uses the diagonal of A rather than solving the full
# system, instead of building and inverting the whole Jacobian. Newton_Raphson() still
# uses the full Jacobian (via Deriv()) and inverts it, same as before.

from MyLib import *

def Copy_Mat(A):
	return [row[:] for row in A]

def Solve(A, b):
	A = Copy_Mat(A)
	n = len(A)
	symmetric = True
	for i in range(n):
		for j in range(n):
			if A[i][j] != A[j][i]:
				symmetric = False
	if symmetric:
		L = Cholesky(A)
		y = forward(L, b)
		Lt = Cholesky_read(Cholesky_compress(L), 1)
		x = backward(Lt, y)
	else:
		LUmat, perm = LU(A)
		bperm = [b[p] for p in perm]
		y = forward(LUmat, bperm)
		x = backward(LUmat, y)
	return x

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

def GaussPRNG(seed=42.0, iter=1):		# Box-Muller, built on PRNG's uniform stream
	u = list(PRNG(seed, 2*iter))
	for i in range(iter):
		u1 = u[2*i]
		u2 = u[2*i+1]
		yield math.sqrt(-2*math.log(u1)) * math.cos(2*math.pi*u2)

def Deriv(f, x, i=None, h=1e-6):		# Secant-method derivative estimate at point x, step h
										# f: scalar function f(x), or one component F[i] of a multivariable system
										# x: point to evaluate at (float for 1D, list for multivariable)
										# i: index of the variable to perturb (required if x is a list)
	if type(x) is not list:			# 1D case
		return (f(x + h) - f(x)) / h
	else:								# multivariable case: partial derivative w.r.t. x[i]
		x_fwd = x[:]
		x_fwd[i] += h
		return (f(x_fwd) - f(x)) / h

def Newton_Raphson(F, x0, h=1e-6, tolerance=1e-6, max_iter=10000):
	# F is a single function f(x) for 1D, or a list of functions [F1, F2, ...] for multivariable
	# x0 is the initial guess: float for 1D, list of floats for multivariable
	# h is the step used by Deriv() to estimate derivatives/the Jacobian
	iters = 0
	if type(x0) is not list:			# 1D case
		x = x0
		while iters < max_iter:
			fx = F(x)
			dfx = Deriv(F, x, h=h)
			x_new = x - fx/dfx
			iters += 1
			if abs(x_new - x) < tolerance:
				break
			x = x_new
		if iters == max_iter:
			print(f"Newton-Raphson did not converge in {max_iter} iterations")
		else:
			print(f"Newton-Raphson converged in {iters} iterations")
		return x_new
	else:								# multivariable case
		n = len(x0)
		x = x0[:]
		while iters < max_iter:
			fx = [F[i](x) for i in range(n)]
			Jx = [[Deriv(F[i], x, j, h) for j in range(n)] for i in range(n)]		# Jacobian, estimated entry by entry
			Jinv = invert(Jx)
			dx = [sum(Jinv[i][j]*fx[j] for j in range(n)) for i in range(n)]
			x_new = [x[i] - dx[i] for i in range(n)]
			iters += 1
			mag = math.sqrt(sum((x_new[i]-x[i])**2 for i in range(n)))
			x = x_new
			if mag < tolerance:
				break
		if iters == max_iter:
			print(f"Newton-Raphson did not converge in {max_iter} iterations")
		else:
			print(f"Newton-Raphson converged in {iters} iterations")
		return x_new

def Fixed_Point(F, x0, h=1e-6, tolerance=1e-6, max_iter=10000):
	# Same signature/derivative-free spirit as Newton_Raphson() above, but the multivariable
	# case only ever uses each equation's own self-derivative (Jacobi-style), never the full
	# Jacobian or a matrix inversion. In 1D this reduces to the same update as
	# Newton_Raphson(), since a 1x1 Jacobian IS just the self-derivative.
	iters = 0
	if type(x0) is not list:			# 1D case
		x = x0
		while iters < max_iter:
			fx = F(x)
			dfx = Deriv(F, x, h=h)
			x_new = x - fx/dfx
			iters += 1
			if abs(x_new - x) < tolerance:
				break
			x = x_new
		if iters == max_iter:
			print(f"Fixed Point did not converge in {max_iter} iterations")
		else:
			print(f"Fixed Point converged in {iters} iterations")
		return x_new
	else:								# multivariable case, diagonal/Jacobi-style
		n = len(x0)
		x = x0[:]
		while iters < max_iter:
			x_new = x[:]
			for i in range(n):
				fx_i = F[i](x)					# uses OLD x throughout the sweep, like Jacobi() in MyLib.py
				dfx_i = Deriv(F[i], x, i, h)	# self-derivative only: dF_i/dx_i, ignoring cross terms
				x_new[i] = x[i] - fx_i/dfx_i
			diff = [x_new[i]-x[i] for i in range(n)]
			mag = math.sqrt(sum(d**2 for d in diff))
			iters += 1
			x = x_new
			if mag < tolerance:
				break
		if iters == max_iter:
			print(f"Fixed Point did not converge in {max_iter} iterations")
		else:
			print(f"Fixed Point converged in {iters} iterations")
		return x_new
