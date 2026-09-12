from MyLib import *
q = [[1.0, -1.0, -7.0, 1.0, 6.0],	# x^4 - x^3 - 7x^2 + x + 6 
     [1.0, 0.0, -5.0, 0.0, 4.0], 	# x^4 - 5x^2 + 4
     [2.0, 0.0, -19.5, 0.5, 13.5, -4.5]]# 2x^5 - 19.5x^2 + 0.5x^2 + 13.5x - 4.5
for i in q:
	print(f"Roots of | {poly_str(i)} | : {Laguerre(i)}")
