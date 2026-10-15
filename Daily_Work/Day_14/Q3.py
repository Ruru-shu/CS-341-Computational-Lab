from MyLib import *

# Q3: y'' = -2 on [0, 1], y(0) = 1, y(1) = 3, by the finite element method
def f(x):
	return -2.0

def exact(x):
	return -x**2 + 3*x + 1

N = 4
xs, ys = Finite_Elim(f, 0.0, 1.0, 1.0, 3.0, N)

print("x        FEM        exact      error")
for i in range(N+1):
	print(f"{xs[i]:.3f}    {ys[i]:.6f}   {exact(xs[i]):.6f}   {abs(ys[i] - exact(xs[i])):.2e}")

worst = 0.0
for i in range(N):
	mid = (xs[i] + xs[i+1])/2
	err = abs((ys[i] + ys[i+1])/2 - exact(mid))		# straight line of the element against the exact curve
	if err > worst:
		worst = err
print("Largest Error between the Nodes (Element Midpoints) =", worst)

xs2, ys2 = Finite_Elim(f, 0.0, 1.0, 1.0, 3.0, 2)
xe = [i/100 for i in range(101)]
ye = [exact(x) for x in xe]

plt.figure()
plt.plot(xe, ye, "k-", label="analytical  -x^2 + 3x + 1")
plt.plot(xs2, ys2, "bo--", label="FEM, 2 elements")
plt.plot(xs, ys, "rs--", label="FEM, 4 elements")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Q3: Finite Element Solution vs Analytical")
plt.legend()
plt.grid(True)
plt.show()
