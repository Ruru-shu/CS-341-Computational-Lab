from MyLib import *
def f1(x):
	return float(1/x)
def f2(x):
	return (x*math.cos(x))
Midpoint_1 = Midpoint_Int(f1, 1, 2, 289)	# N was determined in all cases using the corresponding Error Bound formula and then used ceil() or rounded to nearest even number
Simpson_1 = Simpson_Int(f1, 1, 2, 20)
Midpoint_2 = Midpoint_Int(f2, 0, math.pi/2, 610)
Simpson_2 = Simpson_Int(f2, 0, math.pi/2, 22)

print(f"The value of integ{1}{2} 1/x dx:\nMidpoint: {Midpoint_1:.6f} ; Error% = {((0.69314718-Midpoint_1)/0.69314718)*100}%")
print(f"Simpson: {Simpson_1:.6f} ; Error% = {((0.69314718-Simpson_1)/0.69314718)*100}%")
print(f"The value of integ{0}{math.pi/2} x cos x dx:\nMidpoint: {Midpoint_2:.6f} ; Error% = {((math.pi/2 - 1 - Midpoint_2)/(math.pi/2 -1))*100}%")
print(f"Simpson: {Simpson_2:.6f} ; Error% = {((math.pi/2 - 1 - Simpson_2)/(math.pi/2 -1))*100}%")
