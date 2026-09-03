from MyLib import *

A=[]
B=[]

with open("A.txt", "r") as file: #coefficient matrix
	for line in file:
		A.append(list(map(float, line.strip().split())))
with open("B.txt", "r") as file: #values
        for line in file:
                B.append(float(line))

A, perm = LU(A)
B_perm = [B[perm[i]] for i in range(len(B))]
y = forward(A, B_perm)
x = backward(A, y)
print("Solution (a1 to a6):", [round(i, 6) for i in x])
