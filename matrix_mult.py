A = []
B = []
C = []
D = []

with open("asgn0_matA", "r") as file:
	for line in file:
		A.append(list(map(float, line.strip().split())))

with open("asgn0_matB", "r") as file:
	for line in file:
		B.append(list(map(float, line.strip().split())))
with open("asgn0_vecC", "r") as file:
        for line in file:
                C.append(float(line))
with open("asgn0_vecD", "r") as file:
        for line in file:
                D.append(float(line))
AB = []
for i in range(len(A)):
	row = []
	for j in range(len(B[0])):
		sum = 0
		for k in range(len(A[0])):
			sum += A[i][k] * B[k][j]
		row.append(round(sum, 3))
	AB.append(row)
BC = []
for i in range(len(B)):
        row = []
        for j in range(1):
                sum = 0
                for k in range(len(B[0])):
                        sum += B[i][k] * C[k]
                row.append(round(sum, 3))
        BC.append(row)
D_T = [D]
DC = []
for i in range(len(D_T)):
        row = []
        for j in range(1):
                sum = 0
                for k in range(len(D_T[0])):
                        sum += D_T[i][k] * C[k]
                row.append(round(sum, 3))
        DC.append(row)

print("AB")
for i in AB:
	print(i)
print()
print("BC")
for i in BC:
	print(i)
print()
print("D . C")
for i in DC:
	print(i)
