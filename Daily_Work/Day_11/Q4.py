from MyLib import *
def f(x):			# full integrand e^-x * x^4; Gaussian_Int multiplies by e^x internally to extract g(x) = x^4
	return math.exp(-x) * x**4

true_val = 4*3*2

result = Gaussian_Int(f, 0, float('inf'), n=4)
err    = abs(result - true_val) / true_val * 100

print(f"Integral of e^-x * x^4 from 0 to inf")
print(f"Analytical value (4!)  : {true_val}")
print(f"Gauss-Laguerre (n=4)  : {result:.9f}    Error: {err:.6f}%   [4 points]")
