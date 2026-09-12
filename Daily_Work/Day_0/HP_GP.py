t0 = 1.25
cd = 1.5
cr = 0.5
HP_t = t0
GP_t = t0
HP_s = t0
GP_s = t0

for _ in range(14):
	HP_t += cd
	GP_t *= cr
	HP_s += float(1/HP_t)
	GP_s += GP_t
print(f"GP sum: {GP_s} \nHP sum: {HP_s}")
