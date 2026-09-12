from MyLib import *

def f(x):
	return (3*x + math.sin(x) - math.exp(x))

def df(x):
	return 3 + math.cos(x) - math.exp(x)

root_b  = Bisection(f, -1.5, 1.5)
print(f"Bisection root:     {root_b:.6f},  f(root) = {f(root_b):.2e}")

root_rf = Regula_Falsi(f, -1.5, 1.5)
print(f"Regula Falsi root:  {root_rf:.6f},  f(root) = {f(root_rf):.2e}")

root_nr = Newton_Raphson(f, df, 0.0)
print(f"Newton-Raphson root:{root_nr:.6f},  f(root) = {f(root_nr):.2e}")
