# -------------------------------
# PRNG using Linear Congruential Generator (LCG)
# -------------------------------

seed = 12345

def prng():
    global seed
    a = 1103515245
    c = 12345
    m = 2**31

    seed = (a * seed + c) % m
    return seed

# Returns a random integer from 0 to 9
def random_digit():
    return prng() % 10


# -------------------------------
# Matrix Generation
# -------------------------------

def generate_matrix(rows, cols):
    matrix = []

    for i in range(rows):
        row = []
        for j in range(cols):
            row.append(random_digit())
        matrix.append(row)

    return matrix


# -------------------------------
# Save Matrix to File
# -------------------------------

def save_matrix(filename, matrix):
    with open(filename, "w") as file:
        for row in matrix:
            for value in row:
                file.write(str(value) + " ")
            file.write("\n")


# -------------------------------
# Read Matrix from File
# -------------------------------

def load_matrix(filename):
    matrix = []

    with open(filename, "r") as file:
        for line in file:
            row = list(map(int, line.split()))
            matrix.append(row)

    return matrix


# -------------------------------
# Matrix Multiplication
# -------------------------------

def multiply_matrices(A, B):
    rowsA = len(A)
    colsA = len(A[0])
    rowsB = len(B)
    colsB = len(B[0])

    if colsA != rowsB:
        print("Matrix multiplication not possible.")
        return None

    result = []

    for i in range(rowsA):
        row = []

        for j in range(colsB):
            total = 0

            for k in range(colsA):
                total += A[i][k] * B[k][j]

            row.append(total)

        result.append(row)

    return result


# -------------------------------
# Main Program
# -------------------------------

ROWS = 100
COLS = 100

# Generate matrices
matrix1 = generate_matrix(ROWS, COLS)
matrix2 = generate_matrix(ROWS, COLS)

# Save to files
save_matrix("matrix1.txt", matrix1)
save_matrix("matrix2.txt", matrix2)

# Read back into new matrices
new_matrix1 = load_matrix("matrix1.txt")
new_matrix2 = load_matrix("matrix2.txt")

# Multiply
result = multiply_matrices(new_matrix1, new_matrix2)

# Save result
save_matrix("result_matrix.txt", result)

print("Done!")
print("matrix1.txt created")
print("matrix2.txt created")
print("result_matrix.txt created")
