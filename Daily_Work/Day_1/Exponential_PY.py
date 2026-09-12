import math
import matplotlib.pyplot as plt
from MyLib import PRNG

x = list(PRNG(seed=42.0, iter=5000))
y = [-math.log(i) for i in x if i > 0]

plt.hist(y, bins=100)
plt.xlabel("y")
plt.ylabel("Frequency")
plt.title("Exponential distribution from uniform pRNG")
plt.show()
