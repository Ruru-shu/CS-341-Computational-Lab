from MyLib import *
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

n = 5000                        # number of molecules
steps = 10*n                    # relaxation takes ~n steps, so 10n is plenty
A = [1]*n                       # 1 implies the particle is in left
k = n                           # molecules on left, updated incrementally
Left = []
Right = []
time = []

for t in range(steps):
    time.append(t)
    Left.append(k)
    Right.append(n-k)
    i = int(PRNG()*n)           # pick exactly ONE molecule at random
    if A[i]:                    # left -> right
        A[i] = 0
        k -= 1
    else:                       # right -> left
        A[i] = 1
        k += 1

print(f"final Left = {Left[-1]}, Right = {Right[-1]}")
plt.plot(time, Left, label='Left')
plt.plot(time, Right, label='Right')
plt.xlabel("Time")
plt.ylabel("Particles on each side")
plt.legend()
plt.savefig("ehrenfest.png", dpi=120)
