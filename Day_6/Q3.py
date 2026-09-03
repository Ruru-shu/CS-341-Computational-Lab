import math
from MyLib import *

def f(x):			# The entire function is used to parse the function as necessary, and read as such
	return -x - math.cos(x)

a, b = Bracketing(f, 2, 4, beta=0.5)		# The entire function is passed as an argument.
print(f"Bracket found(a,b): [{a:.6f}, {b:.6f}]")
print(f"f(a) = {f(a):.6f},  f(b) = {f(b):.6f}")
print()
root = Regula_Falsi(f, a, b)				# # The entire function is passed as an argument.
print(f"Regula Falsi root: {root:.6f},  f(root) = {f(root):.2e}")
