from MyLib import *

def g(x):                      # f(x) = x^2 - 2x - 3 = 0
	return math.sqrt(2*x + 3)  # g(x) = sqrt(2x+3)

root = Fixed_Point(g, 4.0)
print(f"Fixed Point root: {root:.6f},  f(root) = {root**2 - 2*root - 3:.2e}")
