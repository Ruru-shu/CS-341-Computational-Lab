from MyLib import *

# Q1: beam fixed at both ends, y'' = -(q/(2*Y*I)) * x * (x - L), y(0) = 0, y(L) = 0
L = 2.0
a = 0.05			# side of the square cross-section
q = 5000.0
Y = 2.0e11
I = a**4/12			# moment of inertia of a square section, a*a^3/12
k = q/(2*Y*I)
h = 0.01

def F(x, U):		# U = [y, y']
	return [U[1], -k*x*(x - L)]

def depth(x):		# y at x, using the converged slope
	xs_t, ys_t = RK4(F, 0.0, [0.0, slope], x, h)
	return ys_t[-1][0]

def tilt(x):		# y' at x, using the converged slope
	xs_t, ys_t = RK4(F, 0.0, [0.0, slope], x, h)
	return ys_t[-1][1]

def miss(x):		# zero where the depression is half of the maximum
	return depth(x) - half

xs, shots = Shooting(F, 0.0, L, 0.0, 0.0, 0.0, -0.03, h)
slope = shots[-1][0][1]
print("Converged slope y'(0) =", slope)

x_max = Bisection(tilt, 0.0, L)			# y' = 0 at the maximum depression
y_max = depth(x_max)
half = y_max/2
x_half1 = Bisection(miss, 0.0, x_max)
x_half2 = Bisection(miss, x_max, L)

print("Maximum Depression =", abs(y_max), "m at x =", x_max, "m")
print("Half-Maximum depression =", abs(half), "m at x =", x_half1, "m and x =", x_half2, "m")

plt.figure()
for i in range(len(shots)):
	ys = [p[0] for p in shots[i]]
	plt.plot(xs, ys, label=f"shot {i+1}, y'(0) = {shots[i][0][1]:.5f}")
plt.plot([0.0, L], [0.0, 0.0], "ko", label="Boundary Points")
plt.plot([x_max], [y_max], "r*", markersize=12, label="Maximum")
plt.plot([x_half1, x_half2], [half, half], "gs", label="Half-Maximum")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Q1: Shooting Trajectories of the Beam")
plt.legend()
plt.grid(True)
plt.show()
