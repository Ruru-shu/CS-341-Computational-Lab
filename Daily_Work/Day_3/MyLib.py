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

def LU(A):   #DoLittle Used
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

	return A
