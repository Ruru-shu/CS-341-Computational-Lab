from MyLib import *

F = [lambda v: 3*v[0]**2 - 2*v[0]*v[1] - 3,             # f1(x,y) = 3x^2 - 2xy - 3 = 0
	lambda v: 3*v[1]**2 - 4*v[0]*v[1]]                  # f2(x,y) = 3y^2 - 4xy = 0

J = [[lambda v: 6*v[0] - 2*v[1], lambda v: -2*v[0]],    # 6x -2y ; -2x
	[lambda v: -4*v[1], lambda v: 6*v[1] - 4*v[0]]]     # -4y  ; 6y - 4x

x0 = [2.0, 3.0]
root = Newton_Raphson(F, J, x0)
print(f"Starting guess: {x0}")
print(f"Root: x = {root[0]:.6f}, y = {root[1]:.6f}")
print(f"Verification: f1 = {F[0](root):.2e}, f2 = {F[1](root):.2e}")
