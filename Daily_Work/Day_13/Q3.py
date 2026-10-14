from MyLib import *
import matplotlib.pyplot as plt

k = 1.0
m = 1.0
mu = 0.15
w2 = k/m
h = 0.01

def f(t, Y):			# Y = [x, v]
	return [Y[1], -mu*Y[1] - w2*Y[0]]

ts, ys = RK4(f, 0.0, [1.0, 0.0], 40.0, h)
x = [Y[0] for Y in ys]
v = [Y[1] for Y in ys]
E = [0.5*m*v[i]**2 + 0.5*k*x[i]**2 for i in range(len(ts))]

plt.figure()
plt.plot(ts, x)
plt.xlabel("t")
plt.ylabel("x")
plt.title("x vs t")
plt.grid()

plt.figure()
plt.plot(ts, v)
plt.xlabel("t")
plt.ylabel("v")
plt.title("v vs t")
plt.grid()

plt.figure()
plt.plot(ts, E)
plt.xlabel("t")
plt.ylabel("E")
plt.title("Total energy vs t")
plt.grid()

plt.show()
