from MyLib import *

A = []

with open("matrix", "r") as file: #coefficient matrix
	for line in file:
		A.append(list(map(float, line.strip().split())))

invert(A)

for i in A:
    print(i, "\n")

