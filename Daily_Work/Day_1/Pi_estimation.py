import math
import matplotlib.pyplot as plt

from MyLib import PRNG

x = list(PRNG(42.0, 5000))
y = list(PRNG(9.0, 5000))

k = 0
count = []

for i, j in zip(x, y):
    if (math.sqrt(i**2 + j**2) <= 1.0):
        k += 1
    count.append(k)

x = list(range(21, 5001))
y = [4*count[i]/(i+1) for i in range(20, 5000)]

plt.scatter(x, y)
plt.axhline(math.pi)
plt.xlabel("Number of points")
plt.ylabel("Estimated value of pi")
plt.show()
