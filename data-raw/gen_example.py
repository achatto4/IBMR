import numpy as np
rng = np.random.default_rng(2025)

# ---- DGP: BMI (exposure) -> CAD (primary), T2D (auxiliary) ----
# We want:
#  * true positive BMI->CAD effect (theta_cad ~ +0.4)
#  * T2D a genuine causal/related outcome (theta_t2d ~ +0.5)
#  * a subset of invalid instruments (pleiotropy) SHARED across CAD & T2D  -> high coheterogeneity
#  * clean strong instruments so both IBPRESSO and IBMODE behave
K = 250                          # instruments
theta_cad = 0.40                 # true BMI->CAD
theta_t2d = 0.55                 # true BMI->T2D
N = 3e5                          # GWAS sample size -> se ~ 1/sqrt(N)
se = 1/np.sqrt(N)

# instrument strengths: all strong (|Z|>=5), realistic BMI effect sizes 0.02-0.06
bx_true = rng.uniform(0.02, 0.06, K) * rng.choice([-1,1], K)

# pleiotropy: 30% invalid; of those, a fraction SHARED across both outcomes (D_ov ~ 0.7)
frac_invalid = 0.30
n_inv = int(frac_invalid*K)
inv = rng.permutation(K)[:n_inv]
n_shared = int(0.7*n_inv)
shared = inv[:n_shared]; cad_only = inv[n_shared: n_shared + (n_inv-n_shared)//2]; t2d_only = inv[n_shared + (n_inv-n_shared)//2:]

alpha_cad = np.zeros(K); alpha_t2d = np.zeros(K)
# shared pleiotropy: correlated direct effects on both CAD and T2D (drives coheterogeneity)
sh = rng.normal(0, 0.010, n_shared)
alpha_cad[shared] = sh;  alpha_t2d[shared] = sh*0.9 + rng.normal(0,0.002,n_shared)
alpha_cad[cad_only] = rng.normal(0,0.010,len(cad_only))
alpha_t2d[t2d_only] = rng.normal(0,0.010,len(t2d_only))

# observed betas: by = theta*bx + alpha + noise
bxh  = bx_true + rng.normal(0, se, K)
bycad = theta_cad*bx_true + alpha_cad + rng.normal(0, se, K)
byt2d = theta_t2d*bx_true + alpha_t2d + rng.normal(0, se, K)
sebx = np.full(K, se); secad = np.full(K, se); set2d = np.full(K, se)

# ---- sanity: IVW BMI->CAD should be POSITIVE (near 0.4, inflated by pleiotropy) ----
def ivw(bx,by,sy):
    w=1/sy**2; return np.sum(w*bx*by)/np.sum(w*bx**2)
print(f"IVW BMI->CAD (all)   = {ivw(bxh,bycad,secad):+.3f}   [true {theta_cad}, expect positive, biased up by pleiotropy]")
print(f"IVW BMI->CAD (valid) = {ivw(bxh[~np.isin(np.arange(K),inv)],bycad[~np.isin(np.arange(K),inv)],secad[~np.isin(np.arange(K),inv)]):+.3f}   [should be near {theta_cad}]")
print(f"IVW BMI->T2D (all)   = {ivw(bxh,byt2d,set2d):+.3f}")
print(f"instrument |Z| range = [{np.min(np.abs(bxh/sebx)):.1f}, {np.max(np.abs(bxh/sebx)):.1f}]  (all strong)")
print(f"invalid: {n_inv} ({n_shared} shared) of {K}  -> expect high coheterogeneity, ~{n_inv} outliers")

# save CSV for R packaging
import csv
snps=[f"rs{1000+i}" for i in range(K)]
w=csv.writer(open("sim_bmi_cad_t2d.csv","w"))  # run from data-raw/
w.writerow(["SNP","BetaXG","seBetaXG","Beta_CAD","se_CAD","Beta_T2D","se_T2D"])
for i in range(K):
    w.writerow([snps[i], round(bxh[i],6), round(sebx[i],6), round(bycad[i],6), round(secad[i],6), round(byt2d[i],6), round(set2d[i],6)])
print("wrote sim_bmi_cad_t2d.csv")
