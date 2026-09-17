from MyLib import *
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def f3(x):
    return math.sin(x)**2

# Exact analytical value: integral of sin^2(x) from -1 to 1
# = [x/2 - sin(2x)/4]_{-1}^{1} = 1 - sin(2)/2
exact = 1 - math.sin(2) / 2

a, b = -1, 1

N_values = list(range(10, 101)) + list(range(101, 10001, 10))

FN_values     = []   # Monte Carlo estimate F_N for each N
sigma_f_vals  = []   # sigma_f: sample std dev of f(X_i) values
                     # this is a PROPERTY OF f itself and stabilises quickly
                     # It does NOT decrease with N; it converges to a constant.
sigma_FN_vals = []   # sigma_{F_N}: standard error of the ESTIMATOR F_N = (b-a) * sigma_f / sqrt(N)
                     # this is what DOES decrease like 1/sqrt(N).
                     #      The "decrease in sigma_f ~ 1/sqrt(N)"
                     #      refers to this quantity: how uncertain our estimate F_N
                     #      is, not sigma_f of the samples which is a fixed constant.
for N in N_values:
    samples = [PRNG() for _ in range(N)]
    Xi = [a + (b - a) * xi for xi in samples]
    fi = [f3(x) for x in Xi]

    FN = (b - a) * sum(fi) / N
    mean_f2 = sum(v**2 for v in fi) / N
    mean_f  = sum(fi) / N
    var     = mean_f2 - mean_f**2
    sigma_f = math.sqrt(abs(var))

    sigma_FN = (b - a) * sigma_f / math.sqrt(N)

    FN_values.append(FN)
    sigma_f_vals.append(sigma_f)
    sigma_FN_vals.append(sigma_FN)

print(f"Exact value            : {exact:.6f}")
print(f"Monte Carlo (N=10000)  : {FN_values[-1]:.6f}")

for N, FN in zip(N_values, FN_values):
    if abs(FN - exact) < 5e-5:
        print(f"First N accurate to 4dp: N = {N}, FN = {FN:.6f}")
        break

# ── Plot 1 : F_N vs N ────────────────────────────────────────────────────────
fig1, ax1 = plt.subplots(figsize=(9, 5))

ylim_lo, ylim_hi = exact - 0.20, exact + 0.20
FN_masked = [v if ylim_lo <= v <= ylim_hi else None for v in FN_values]

ax1.plot(N_values, FN_masked,
         color="#2E86AB", lw=0.9, alpha=0.85,
         label=r"$\mathcal{F}_N$ (Monte Carlo estimate)")
ax1.axhline(exact, color="#E84855", lw=1.4, ls="--",
            label=f"Exact = {exact:.6f}")

ax1.set_xlabel("N  (number of samples)", fontsize=12)
ax1.set_ylabel(r"$\mathcal{F}_N$", fontsize=13)
ax1.set_title(r"Monte Carlo integration of $\int_{-1}^{1} \sin^2 x\; dx$", fontsize=13)
ax1.legend(fontsize=10)
ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

ax1.set_ylim(ylim_lo, ylim_hi)
ax1.grid(True, alpha=0.3, lw=0.6)
fig1.tight_layout()
fig1.savefig("Q3_FN_vs_N.png", dpi=150)
plt.show()
print("Plot 1 saved: Q3_FN_vs_N.png")

fig2, ax2 = plt.subplots(figsize=(9, 5))

ax2.plot(N_values, sigma_f_vals,
         color="#F4A261", lw=1.0, alpha=0.90,
         label=r"$\sigma_f$  (std dev of sample $f(X_i)$  —  stabilises to a constant)")

ax2.plot(N_values, sigma_FN_vals,
         color="#2E86AB", lw=1.1, alpha=0.85,
         label=r"$\sigma_{\mathcal{F}_N} = (b{-}a)\,\sigma_f/\sqrt{N}$  (std error of estimator  —  falls as $1/\sqrt{N}$)")


sigma_f_stable = sigma_f_vals[-1]
scale = sigma_f_stable * math.sqrt(N_values[0]) * (b - a)
ref   = [scale / math.sqrt(N) for N in N_values]

ax2.plot(N_values, ref,
         color="#264653", lw=1.2, ls=":", alpha=0.65,
         label=r"$\propto 1/\sqrt{N}$ reference")

ax2.set_xlabel("N  (number of samples)", fontsize=12)
ax2.set_ylabel(r"$\sigma$", fontsize=13)
ax2.set_title(r"$\sigma_f$ (Standard Deviation) and $\sigma_{\mathcal{F}_N}$ (Standard Error) vs $N$",
              fontsize=12)
ax2.legend(fontsize=9, loc="upper right")
ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
ax2.set_ylim(0, sigma_f_stable * 1.25)
ax2.grid(True, alpha=0.3, lw=0.6)
fig2.tight_layout()
fig2.savefig("Q3_sigma_vs_N.png", dpi=150)
plt.show()
print("Plot 2 saved: Q3_sigma_vs_N.png")
