"""Referee: global minimum of the discrete K=I closed-chain energy E = sum_j |q_{j+1}-q_j|^2/(2h) over ALL sections of
the admissible self-motion arc bundle through the rest-posture lift (nodes restricted to a grid in fibre arclength),
by min-plus dynamic programming around the ring. Returns the landscape E*(t0) = min energy of closed chains through
the fibre point t0 at node 0, its local minima, and the optimal chain; optionally with a joint-margin constraint."""
import sys, time
import numpy as np
import rk
m, d = int(sys.argv[1]), float(sys.argv[2])
margin = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
W = int(sys.argv[4]) if len(sys.argv) > 4 else 24
arcs = np.load(f"arcs_m{m}_d{d}.npy", allow_pickle=True)
h = rk.L/m
lo = min(-a["k0"] for a in arcs); hi = max(len(a["Q"])-1-a["k0"] for a in arcs); G = hi-lo+1
P = np.full((m, G, 7), np.nan)
for j, a in enumerate(arcs):
    st = -a["k0"]-lo; P[j, st:st+len(a["Q"])] = a["Q"]
ok = np.isfinite(P[..., 0])
mg = np.where(ok, np.minimum(P-rk.QMIN, rk.QMAX-P).min(axis=2), -1)
ok &= mg >= margin
t0 = time.time()
S = np.full((G, G), np.inf); S[np.arange(G), np.arange(G)] = np.where(ok[0], 0.0, np.inf)
back = np.zeros((m, G, G), dtype=np.int8)
for j in range(m):
    A, B = P[j], P[(j+1) % m]; okA, okB = ok[j], ok[(j+1) % m]
    Snew = np.full((G, G), np.inf); arg = np.zeros((G, G), dtype=np.int8)
    for o in range(-W, W+1):
        ka = np.arange(max(0, -o), min(G, G-o)); kb = ka+o
        c = np.sum((B[kb]-A[ka])**2, 1)/(2*h)
        c = np.where(okA[ka] & okB[kb] & np.isfinite(c), c, np.inf)
        cand = S[:, ka] + c[None, :]
        cur = Snew[:, kb]; better = cand < cur
        cur[better] = cand[better]; Snew[:, kb] = cur
        ag = arg[:, kb]; ag[better] = o; arg[:, kb] = ag
    back[j] = arg; S = Snew
prof = S[np.arange(G), np.arange(G)]
i0 = int(np.argmin(prof))
ks = [i0]; k = i0
for j in range(m-1, -1, -1):
    k = k-int(back[j][i0, k]); ks.append(k)
ks = ks[::-1]; assert ks[0] == i0, (ks[0], i0)
Qopt = np.array([P[j, ks[j]] for j in range(m)])
E = np.sum((np.roll(Qopt, -1, 0)-Qopt)**2)/(2*h)
tt = (np.arange(G)+lo)*d
mgopt = np.minimum(Qopt-rk.QMIN, rk.QMAX-Qopt).min()
print(f"m={m} d={d} margin>={margin} W={W}: DP global min E = {prof[i0]:.4f} (recomputed {E:.4f}) at t0 = {tt[i0]:.2f}; "
      f"chain t range [{(np.array(ks[:-1])+lo).min()*d:.2f},{(np.array(ks[:-1])+lo).max()*d:.2f}]; "
      f"max|offset| {np.max(np.abs(np.diff(ks)))} of W={W}; min joint margin {mgopt:.4f}; {time.time()-t0:.0f}s")
fin = np.isfinite(prof)
loc = [i for i in range(1, G-1) if fin[i-1] and fin[i] and fin[i+1] and prof[i] <= prof[i-1] and prof[i] <= prof[i+1]]
print("local minima of E*(t0) (t0 = fibre arclength from the rest-lift node at s=0, E*):",
      [(round(tt[i], 2), round(prof[i], 4)) for i in loc])
np.savez(f"dp2_m{m}_d{d}_mg{margin}.npz", prof=prof, tt=tt, Qopt=Qopt, ks=np.array(ks))
