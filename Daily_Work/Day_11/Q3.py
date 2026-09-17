from MyLib import *
def f(x):			# full integrand; Gaussian_Int multiplies by e^x internally to extract g(x) = 1/(1+x)
	return math.exp(-x) / (1 + x)

true_val = 0.59635

result = Gaussian_Int(f, 0, float('inf'), n=5)
err    = abs(result - true_val) / true_val * 100

print(f"Integral of e^-x/(1+x) from 0 to inf")
print(f"True value             : {true_val}")
print(f"Gauss-Laguerre (n=5)  : {result:.9f}    Error: {err:.6f}%   [5 points]")
