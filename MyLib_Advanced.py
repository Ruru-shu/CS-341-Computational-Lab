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
