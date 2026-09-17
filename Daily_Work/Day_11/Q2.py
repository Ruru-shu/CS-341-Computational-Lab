from MyLib import *

def f(x):
	return x**2 / (1 + x**4)

true_val = 0.487495494

gauss_result  = Gaussian_Int(f, -1, 1, 4)
simp_result   = Simpson_Int(f, -1, 1, 60)

gauss_err = abs(gauss_result - true_val) / true_val * 100
simp_err  = abs(simp_result  - true_val) / true_val * 100

print(f"Integral of x^2/(1+x^4) from -1 to 1")
print(f"True value                : {true_val}")
print(f"Gauss-Legendre (n=4)      : {gauss_result:.9f}    Error: {gauss_err:.6f}%   [4 points]")
print(f"Simpson's 1/3-rule (N=60) : {simp_result:.9f}    Error: {simp_err:.6f}%   [60 points]")
