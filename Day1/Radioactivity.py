import math
import matplotlib.pyplot as plt
from MyLib import PRNG

lambda_A = 0.0347
lambda_B = 0.0198

NA = 500
NB = 0
NC = 0

dt = 1.0
steps = 300

rand_A = list(PRNG(seed=42.0, iter=steps*600))
rand_B = list(PRNG(seed=17.0, iter=steps*600))

time = [0]
NA_list = [NA]
NB_list = [NB]
NC_list = [NC]

idx = 0
for step in range(steps):
    dA = 0
    for _ in range(NA):
        if rand_A[idx] < lambda_A * dt:
            dA += 1
        idx += 1

    dB = 0
    for _ in range(NB):
        if rand_B[idx % len(rand_B)] < lambda_B * dt:
            dB += 1
        idx += 1

    NA -= dA
    NB += dA - dB
    NC += dB

    time.append((step+1)*dt)
    NA_list.append(NA)
    NB_list.append(NB)
    NC_list.append(NC)

plt.plot(time, NA_list, label='NA')
plt.plot(time, NB_list, label='NB')
plt.plot(time, NC_list, label='NC')
plt.xlabel("Time")
plt.ylabel("Number of nuclei")
plt.legend()
plt.show()
