from MyLib import *
A=[]
B=[]
with open("A2.txt", "r") as file: #coefficients
	for line in file:
		A.append(list(map(float, line.strip().split())))
with open("B2.txt", "r") as file: #values
        for line in file:
                B.append(float(line))

x_j = Jacobi([row[:] for row in A], B, tolerance=1e-6)
x_j = [round(v, 6) for v in x_j]
print("Solved via Jacobi\nx:", x_j)
print()
x_gs = Gauss_Seidel([row[:] for row in A], B, tolerance=1e-6)
x_gs = [round(v, 6) for v in x_gs]
print("Solved via Gauss_Seidel\nx:", x_gs)
