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
def PRNG(seed=42.0, iter=1):                         #Linear Congruential Generator; Seed and number of iterations are passed in function header
	a = 1103515245.0
	c = 12345.0
	m = 32768.0
	for i in range(iter):
		seed = (a*seed +c) % m
		yield float(seed/m)
