from MyLib import *
A=[]
B=[]
with open("GJE_Q.txt", "r") as file: #coefficients
	for line in file:
		A.append(list(map(float, line.strip().split())))
with open("GJE_S.txt", "r") as file: #values
        for line in file:
                B.append(float(line))

x = GJE(A, B)
print("x, y, z =", x)
