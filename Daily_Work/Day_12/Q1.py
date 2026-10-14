from MyLib import *
import matplotlib.pyplot as plt

def f1(x, y):
	return y - x**2

def exact1(x):
	return x**2 + 2*x + 2 - 2*math.exp(x)

def f2(x, y):
	return (x + y)**2

def exact2(x):
	return math.tan(x + math.pi/4) - x

steps = [0.02, 0.1, 0.3, 0.5]

# ODE 1 : x from 0 to 2
xe = [2*i/100 for i in range(101)]
ye = [exact1(x) for x in xe]
for h in steps:
	plt.figure()
	plt.plot(xe, ye, "k", linewidth=2, label="Analytical")
	x, y = Euler_Fwd(f1, 0, 0, 2, h)
	plt.plot(x, y, ".-", label="Forward Euler")
	x, y = Euler_Bkwd(f1, 0, 0, 2, h)
	plt.plot(x, y, ".-", label="Backward Euler")
	x, y = Predictor_Corrector(f1, 0, 0, 2, h)
	plt.plot(x, y, ".-", label="Predictor-Corrector")
	plt.title(f"ODE 1: dy/dx = y - x^2,  step h = {h}")
	plt.xlabel("x")
	plt.ylabel("y")
	plt.legend()

# ODE 2 : x from 0 to pi/5
xe = [(math.pi/5)*i/100 for i in range(101)]
ye = [exact2(x) for x in xe]
for h in steps:
	plt.figure()
	plt.plot(xe, ye, "k", linewidth=2, label="Analytical")
	x, y = Euler_Fwd(f2, 0, 1, math.pi/5, h)
	plt.plot(x, y, ".-", label="Forward Euler")
	x, y = Euler_Bkwd(f2, 0, 1, math.pi/5, h)
	plt.plot(x, y, ".-", label="Backward Euler")
	x, y = Predictor_Corrector(f2, 0, 1, math.pi/5, h)
	plt.plot(x, y, ".-", label="Predictor-Corrector")
	plt.title(f"ODE 2: dy/dx = (x + y)^2,  step h = {h}")
	plt.xlabel("x")
	plt.ylabel("y")
	plt.legend()

plt.show()
