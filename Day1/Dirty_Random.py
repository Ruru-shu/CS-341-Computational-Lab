import matplotlib.pyplot as plt
def Dirty_Random(seed=0.1, iter=1, c=-0.13):
	x0 = seed
	for i in range(iter):
		x0 = c*x0*(1-x0)
		yield (x0)
k = 5
c_vals = [3.9, 3.7, 3.8]

for c in c_vals:
	x = list(Dirty_Random(seed=0.1, iter=1000, c=c))
	y = list(Dirty_Random(seed=0.1, iter=1005, c=c))[5:]

	plt.scatter(y, x)
	plt.xlabel("x[i+k]")
	plt.ylabel("x[i]")
	plt.title(f"c = {c}, k = {k}")
	plt.show()
