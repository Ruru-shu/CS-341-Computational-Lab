k = 0
m = 2
sum = 0
fact = 1
n = 8
while True:
	if m%2 == 0:
		sum +=m
		m+=2
		k+=1
	if k==20:
		break
while (n!=1):
	fact*=n
	n-=1
print(f"Sum of first 20 even numbers: {sum}")	
print(f"Factorial of 8 is: {fact}")
