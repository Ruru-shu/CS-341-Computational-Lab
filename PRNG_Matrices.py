seed = 42.0

def PRNG():
    global seed
    a = 1664525.0
    c = 1013904223.0
    m = 2**32
    seed = (a * seed + c) % m
    return float(seed/m)

A = []
B = []
for i in range(100):
	row = []
	row2 = []
	for j in range(100):
		row.append(int(PRNG()*10))
		row2.append(int(PRNG()*10))
	A.append(row)
	B.append(row2)

with open("matrix_A.txt", "w") as file:
	for i in range(100):
		for j in range(100):
			file.write(f"{A[i][j]} ")
		file.write("\n")

with open("matrix_B.txt", "w") as file:
	for i in range(100):
		for j in range(100):
			file.write(f"{B[i][j]} ")
		file.write("\n")

C = []
D = []

with open("matrix_A.txt", "r") as file:
	for line in file:
		C.append(list(map(int, line.strip().split())))

with open("matrix_B.txt", "r") as file:
	for line in file:
		D.append(list(map(int, line.strip().split())))

E = []
for i in range(len(C)):
	crow = []
	for j in range(len(C)):
		row = 0
		for k in range(len(C)):
			row += C[i][k] * D[k][j]
		crow.append(row)
	E.append(crow)

for i in range(len(E)):
	for j in range(len(E[0])):
		print(f"{E[i][j]} ", end='')
	print()
