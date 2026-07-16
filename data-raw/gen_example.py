# Generate the simulated illustrative dataset for IBMR (seed 2025).
# One exposure, one primary outcome, and TWO candidate auxiliary outcomes:
#   Auxiliary_1 shares most of its invalid-instrument (pleiotropic) structure
#   with the primary outcome (high coheterogeneity, the informative auxiliary),
#   Auxiliary_2 shares little (low coheterogeneity).
# Instruments are all strong so both IBMODE and IBPRESSO behave stably.
import numpy as np, csv
rng = np.random.default_rng(2025)

K = 250
theta_primary = 0.40          # true exposure -> primary-outcome effect
N = 3e5; se = 1/np.sqrt(N)
bx_true = rng.uniform(0.02, 0.06, K) * rng.choice([-1, 1], K)

n_inv = int(0.30*K); inv = rng.permutation(K)[:n_inv]

a_primary = np.zeros(K); a_aux1 = np.zeros(K); a_aux2 = np.zeros(K)
# Auxiliary_1: 70% of the primary's invalid set is SHARED -> high coheterogeneity
n_sh1 = int(0.70*n_inv); shared1 = inv[:n_sh1]
sh = rng.normal(0, 0.010, n_sh1)
a_primary[shared1] = sh
a_aux1[shared1]    = 0.9*sh + rng.normal(0, 0.002, n_sh1)
# primary-only invalid
prim_only = inv[n_sh1:]
a_primary[prim_only] = rng.normal(0, 0.010, len(prim_only))
# Auxiliary_1 also has its own independent invalid instruments
a_inv1 = rng.permutation(K)[:n_inv]; a_aux1[a_inv1] += rng.normal(0, 0.004, n_inv)
# Auxiliary_2: invalid instruments mostly DISJOINT from primary -> low coheterogeneity
a_inv2 = rng.permutation(K)[:n_inv]; a_aux2[a_inv2] = rng.normal(0, 0.010, n_inv)

bxh   = bx_true + rng.normal(0, se, K)
byp   = theta_primary*bx_true + a_primary + rng.normal(0, se, K)
bya1  = 0.55*bx_true + a_aux1 + rng.normal(0, se, K)   # aux1 also causally related
bya2  = 0.30*bx_true + a_aux2 + rng.normal(0, se, K)
sebx = np.full(K, se); sep = np.full(K, se); sa1 = np.full(K, se); sa2 = np.full(K, se)

# quick coheterogeneity proxy
def cohet(bx,by1,by2,sy1,sy2):
    th1=by1/bx; th2=by2/bx; s1=(sy1/bx)**2; s2=(sy2/bx)**2
    w=1/np.sqrt(s1*s2); w/=w.sum()
    d1=th1-np.sum(w*th1); d2=th2-np.sum(w*th2)
    return np.sum(w*(d1*d2))/np.sqrt(np.sum(w*d1**2)*np.sum(w*d2**2))
print(f"coheterogeneity(primary, Auxiliary_1) = {cohet(bxh,byp,bya1,sep,sa1):.2f}  [expect HIGH]")
print(f"coheterogeneity(primary, Auxiliary_2) = {cohet(bxh,byp,bya2,sep,sa2):.2f}  [expect LOW]")

w=csv.writer(open("sim_ibmr_example.csv","w"))
w.writerow(["SNP","BetaXG","seBetaXG","Beta_Primary","se_Primary","Beta_Aux1","se_Aux1","Beta_Aux2","se_Aux2"])
for i in range(K):
    w.writerow([f"rs{1000+i}", round(bxh[i],6),round(sebx[i],6), round(byp[i],6),round(sep[i],6),
                round(bya1[i],6),round(sa1[i],6), round(bya2[i],6),round(sa2[i],6)])
print("wrote sim_ibmr_example.csv")
