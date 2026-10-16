"""Referee: DP global minimum (as dp2.py) of barrier-modified energies on the main admissible arc bundle:
 mode 'metric': conformal barrier metric K = phi^2 I, phi = 1 + c sum_i 1/((q_i-q_i^-)(q_i^+-q_i)) (theory Prop. 8.1, c_Sigma = 0),
                discrete E = sum_j phi(q_{j+1/2})^2 |q_{j+1}-q_j|^2/(2h)
 mode 'pot'   : barrier potential B = c sum_i [(q_i-q_i^-)^-2 + (q_i^+-q_i)^-2] (theory Prop. 8.2), E + h sum_j B(q_j).
Prints the global minimiser's location (t0 at node 0) and the local minima of the landscape E*(t0)."""
import sys, time
import numpy as np, rk
mode, c = sys.argv[1], float(sys.argv[2]); m, d, W = 50, 0.01, 24
arcs = np.load(f"arcs_m{m}_d{d}.npy", allow_pickle=True); h = rk.L/m
lo = min(-a["k0"] for a in arcs); hi = max(len(a["Q"])-1-a["k0"] for a in arcs); G = hi-lo+1
P = np.full((m, G, 7), np.nan)
for j, a in enumerate(arcs):
    st = -a["k0"]-lo; P[j, st:st+len(a["Q"])] = a["Q"]
ok = np.isfinite(P[..., 0])
def phi(X): return 1 + c*np.sum(1/((X-rk.QMIN)*(rk.QMAX-X)), axis=-1)
def Bpot(X): return c*np.sum((X-rk.QMIN)**-2 + (rk.QMAX-X)**-2, axis=-1)
nodecost = np.where(ok, h*Bpot(np.nan_to_num(P, nan=0.0)) if mode == "pot" else 0.0, np.inf)
S = np.full((G, G), np.inf); S[np.arange(G), np.arange(G)] = nodecost[0]
for j in range(m):
    A, B = P[j], P[(j+1) % m]; okA, okB = ok[j], ok[(j+1) % m]; nc = nodecost[(j+1) % m]
    Snew = np.full((G, G), np.inf)
    for o in range(-W, W+1):
        ka = np.arange(max(0, -o), min(G, G-o)); kb = ka+o
        dq = B[kb]-A[ka]; w = phi(0.5*(A[ka]+B[kb]))**2 if mode == "metric" else 1.0
        cc = w*np.sum(dq**2, 1)/(2*h)
        cc = np.where(okA[ka] & okB[kb] & np.isfinite(cc), cc, np.inf)
        if j < m-1: cc = cc + nc[kb]
        cand = S[:, ka]+cc[None, :]; Snew[:, kb] = np.minimum(Snew[:, kb], cand)
    S = Snew
prof = S[np.arange(G), np.arange(G)]; tt = (np.arange(G)+lo)*d; i0 = int(np.nanargmin(prof))
fin = np.isfinite(prof)
loc = [i for i in range(1, G-1) if fin[i-1] and fin[i] and fin[i+1] and prof[i] <= prof[i-1] and prof[i] <= prof[i+1]]
print(f"{mode} c={c}: global min {prof[i0]:.4f} at t0 = {tt[i0]:.2f}; landscape minima:", [(round(tt[i], 2), round(prof[i], 4)) for i in loc], flush=True)
