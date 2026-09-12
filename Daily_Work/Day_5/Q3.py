from MyLib import *
A=[]
B=[]
with open("A3.txt", "r") as file: #coefficients
	for line in file:
		A.append(list(map(float, line.strip().split())))
with open("B3.txt", "r") as file: #values
        for line in file:
                B.append(float(line))

x_gs = Gauss_Seidel([row[:] for row in A], B)
x_gs = [round(v, 6) for v in x_gs]
print("Solved via Gauss_Seidel:\n", x_gs)
print()
x_sor = Gauss_Seidel([row[:] for row in A], B, 1.57)
x_sor = [round(v, 6) for v in x_sor]
print("Solved via SOR(w=1.57):", x_sor)
