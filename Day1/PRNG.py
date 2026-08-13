import matplotlib.pyplot as plt
def PRNG(seed=42.0, iter=1):                         #Linear Congruential Generator; Seed and number of iterations are passed in function header
	a = 1103515245.0
	c = 12345.0
	m = 32768.0
	for i in range(iter):
		seed = (a*seed +c) % m
		yield float(seed/m)

x = list(PRNG(42.0, 1000))
y = list(PRNG(42.0, 1005))[5:] 

plt.scatter(y, x)
plt.xlabel("x")
plt.ylabel("y")
plt.show()
