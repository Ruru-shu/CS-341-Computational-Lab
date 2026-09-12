import math
from MyLib import *

def f(x):			# The entire function is used to parse the function as necessary, and read as such
	return math.log(x/2) - math.sin(5*x/2)

root_bis = Bisection(f, 1.5, 3.0)			# The entire function is passed as an argument.
print(f"Bisection root:    {root_bis:.6f},  f(root) = {f(root_bis):.2e}")
root_rf = Regula_Falsi(f, 1.5, 3.0)			# # The entire function is passed as an argument.
print(f"Regula Falsi root: {root_rf:.6f},  f(root) = {f(root_rf):.2e}")
