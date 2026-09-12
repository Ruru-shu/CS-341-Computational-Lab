from MyLib import *

def f1(x):
	return float(1/x)

def f2(x):
	return float(x*math.cos(x))

def f3(x):
	return x*math.atan(x)

exact1 = 0.69314718
exact2 = math.pi/2 - 1
exact3 = math.pi/4 - 0.5

N_values = [4, 8, 15, 20]

def print_table(label, f, a, b, exact):
	print(f"\nFor {label}  ;  Exact = {exact:.8f}")
	print(f"{'N':<6} {'Trap Value':<18} {'Trap Error%':<18} {'Mid Value':<18} {'Mid Error%':<18}")
	print("-" * 78)
	for N in N_values:
		trap = Trapezoidal_Int(f, a, b, N)
		mid  = Midpoint_Int(f, a, b, N)
		trap_err = abs(trap - exact)/abs(exact)*100
		mid_err  = abs(mid  - exact)/abs(exact)*100
		print(f"{N:<6} {trap:<18.8f} {trap_err:<18.6f} {mid:<18.8f} {mid_err:<18.6f}")

print_table(r"\int_{1}^{2} (1/x) dx",f1, 1, 2, exact1)
print_table(r"\int_{0}^{pi/2} x*cos(x) dx", f2, 0, math.pi/2, exact2)
print_table(r"\int_{0}^{1} x*arctan(x) dx", f3, 0, 1, exact3)
