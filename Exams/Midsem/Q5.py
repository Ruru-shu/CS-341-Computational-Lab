from MyLib import *

def f1(u):
	return math.exp(-u) * (u**3) / 4

r_avg = Gaussian_Int(f1, 0, float('inf'), n=5)
print(f"Q1: <r> = {r_avg:.4f} a0")

