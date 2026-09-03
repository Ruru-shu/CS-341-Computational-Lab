n = 10
A = [[0.0]*n for i in range(n)]

for i in range(n):
	A[i][i] = 2.0
	if i > 0:
		A[i][i-1] = -1.0
	if i < n-1:
		A[i][i+1] = -1.0

with open("A3.txt", "w") as f:
	for row in A:
		line = " ".join(str(val) for val in row)
		f.write(line + "\n")

with open("B3.txt", "w") as f:
	for i in range(n):
		f.write("1.0\n")
