from MyLib import *
import matplotlib.pyplot as plt

a = 0.2
b = 0.2
c = 5.7

def f(t, Y):			# Y = [x, y, z]
	return [-Y[1] - Y[2], Y[0] + a*Y[1], b + Y[2]*(Y[0] - c)]

ts, ys = RK4(f, 0.0, [1.0, 0.0, 0.0], 100.0, 0.01)
x = [Y[0] for Y in ys]
y = [Y[1] for Y in ys]
z = [Y[2] for Y in ys]

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
ax.plot(x, y, z, linewidth=0.5)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("Rossler system")
plt.show()
