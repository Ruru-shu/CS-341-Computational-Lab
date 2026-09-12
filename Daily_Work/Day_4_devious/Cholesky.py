from MyLib import *

A=[]
B=[]
with open("3A.txt", "r") as file: #coefficient matrix
	for line in file:
		A.append(list(map(float, line.strip().split())))
with open("3B.txt", "r") as file: #values
        for line in file:
                B.append(float(line))

A = Cholesky(A)
n = len(A)

LT = [[0.0]*n for i in range(n)]
for i in range(n):
	for j in range(n):
		LT[i][j] = A[j][i]

y = forward(A, B, unit_diag=False)
x = backward(LT, y)
print("Solution (x1 to x4):", [round(i, 6) for i in x])
