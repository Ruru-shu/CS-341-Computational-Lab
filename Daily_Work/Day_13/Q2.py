from MyLib import *
import matplotlib.pyplot as plt
def f(x, y):
	return (x + y)**2

def exact(x):
	return math.tan(x + math.pi/4) - x

x_end = math.pi/5

xe = [x_end*i/500 for i in range(501)]
ye = [exact(x) for x in xe]
plt.plot(xe, ye, "k-", label="Analytical")

markers = ["o", "s"]
for h, m in zip([0.05, 0.15], markers):
	xs, ys = RK4(f, 0.0, 1.0, x_end, h)
	err = max(abs(ys[i] - exact(xs[i])) for i in range(len(xs)))
	print(f"h = {h}: max |RK4 - analytical| = {err:.3e}")
	plt.plot(xs, ys, m, label=f"RK4, h = {h}")

plt.xlabel("x")
plt.ylabel("y")
plt.title("dy/dx = (x + y)^2, y(0) = 1")
plt.legend()
plt.grid()
plt.show()
