from MyLib import *

def f(x):
	return (x-5)*math.exp(x) + 5

def df(x):
	return (x-5)*math.exp(x) + math.exp(x)

root_nr = Newton_Raphson(f, df, 5.0)
print(f"Newton-Raphson root:{root_nr},  f(root) = {f(root_nr)}")

b = (6.626*3*0.001)/(1.381*root_nr)
print(f"The value of Wein's constant(b): {b:.4f}")
