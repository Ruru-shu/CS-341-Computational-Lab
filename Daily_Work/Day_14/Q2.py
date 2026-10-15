from MyLib import *

# Q2: heat conduction in a thin rod, T'' + alpha*(Ta - T) = 0, T(0) = 40, T(L) = 200
L = 10.0
alpha = 0.01
Ta = 20.0
h = 0.05

def F(x, U):		# U = [T, T']
	return [U[1], alpha*(U[0] - Ta)]

def temp(x):		# T at x, using the converged slope
	xs_t, ys_t = RK4(F, 0.0, [40.0, slope], x, h)
	return ys_t[-1][0]

def miss(x):		# zero where T = 100
	return temp(x) - 100.0

xs, shots = Shooting(F, 0.0, L, 40.0, 200.0, 0.0, 30.0, h)
slope = shots[-1][0][1]
print("Converged slope T'(0) =", slope)

x100 = Bisection(miss, 0.0, L)
print("T = 100 C at x =", x100, "m")

plt.figure()
for i in range(len(shots)):
	ys = [p[0] for p in shots[i]]
	plt.plot(xs, ys, label=f"shot {i+1}, T'(0) = {shots[i][0][1]:.5f}")
plt.plot([0.0, L], [40.0, 200.0], "ko", label="boundary points")
plt.plot([x100], [100.0], "r*", markersize=12, label="T = 100 C")
plt.xlabel("x (m)")
plt.ylabel("T (C)")
plt.title("Q2: Shooting Trajectories of the Rod Temperature")
plt.legend()
plt.grid(True)
plt.show()
