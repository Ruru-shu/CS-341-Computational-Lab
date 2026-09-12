import copy
from MyLib import *

A = []
with open("A.txt", "r") as file:
	for line in file:
		A.append(list(map(float, line.strip().split())))

A_copy = copy.deepcopy(A)
A = LU(A)

n = len(A)
L = [[0.0]*n for i in range(n)]
U = [[0.0]*n for i in range(n)]

for i in range(n):
	for j in range(n):
		if i == j:
			L[i][j] = 1.0
			U[i][j] = A[i][j]
		elif i > j:
			L[i][j] = A[i][j]
		else:
			U[i][j] = A[i][j]

LU = [[sum(L[i][k]*U[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

print("L:")
for i in L:
	print(i)
print()
print("U:")
for i in U:
	print(i)
print()
print("L.U:")
for i in LU:
	print(i)
print()
print("A originally:")
for i in A_copy:
	print(i)

print("The matrices might NOT be exactly equivalent due to row pivoting, so the rows are swapped. But they're identical otherwise.")
