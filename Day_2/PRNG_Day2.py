import matplotlib.pyplot as plt

def PRNG(seed=42.0, iter=1):
	a = 1103515245.0
	c = 12345.0
	m = 32768.0
	if not hasattr(PRNG, "x0"):
		PRNG.x0 = 42.0
	if seed:
		PRNG.x0 = seed
	for i in range(iter):
		x1 = (a * PRNG.x0 + c) % m
		PRNG.x0 = x1
		yield x1 / m

if __name__ == "__main__":
	x = list(PRNG(42.0, 1000))
	y = list(PRNG(42.0, 1005))[5:]

	plt.scatter(y, x)
	plt.xlabel("x")
	plt.ylabel("y")
	plt.show()
