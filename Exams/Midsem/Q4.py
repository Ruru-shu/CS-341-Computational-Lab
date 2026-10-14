from MyLib import *
def lam(x):
	return x**2

def num(x):
	return x * lam(x)

M = Simpson_Int(lam, 0, 2, N=1000)			# total mass
moment = Simpson_Int(num, 0, 2, N=1000)		# first moment
x_cm = moment / M
print(f"Q2: Mass = {M:.4f}, Moment = {moment:.4f}")
print(f"Q2: x_cm = {x_cm:.4f} m")
