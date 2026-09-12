def delta(i, j):
	if (i==j):
		return 1
	else:
		return 0

def perm_matrix(mu, nu):
    matrix = []
    for i in range(3):
        row = []
        for j in range(3):
            element = (delta(i, j)- delta(i, mu) * delta(j, mu)- delta(i, nu) * delta(j, nu)+ delta(i, mu) * delta(j, nu)+ delta(i, nu) * delta(j, mu))
            row.append(element)
        matrix.append(row)
    return matrix

A_xy = perm_matrix(0, 1)
A_zy = perm_matrix(2, 1)
A_xz = perm_matrix(0, 2)

print("A_xy")
for row in A_xy:
    print(row)
print()

print("A_zy")
for row in A_zy:
    print(row)
print()

print("A_xz")
for row in A_xz:
    print(row)
