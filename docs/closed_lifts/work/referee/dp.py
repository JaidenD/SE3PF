"""Global minimum of the discrete closed-chain energy sum_j |q_{j+1}-q_j|^2/(2h) over sections whose nodes lie on
the traced admissible arcs (grid in fibre arclength), by min-plus dynamic programming around the ring.
Also: E*(t0) = min energy of closed chains through fibre point t0 at node 0 (1-D landscape),
and the same with a joint margin constraint."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import rk

here = os.path.dirname(os.path.abspath(__file__))
m = int(sys.argv[1]); d = float(sys.argv[2]); margin = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
W = int(sys.argv[4]) if len(sys.argv) > 4 else 24
arcs = np.load(os.path.join(here, f"arcs_m{m}_d{d}.npy"), allow_pickle=True)
h = rk.L/m
lo = min(-a["k0"] for a in arcs); hi = max(len(a["Q"])-1-a["k0"] for a in arcs)
G = hi - lo + 1
P = np.full((m, G, 7), np.nan)
for j, a in enumerate(arcs):
    Q = a["Q"]; st = -a["k0"]-lo
    P[j, st:st+len(Q)] = Q
mg = np.minimum(P - rk.QMIN, rk.QMAX - P).min(axis=2)
P[~(mg >= margin)] = np.nan
INF = np.inf
S = np.full((G, G), INF); np.fill_diagonal(S, 0.0)
S[np.isnan(P[0, :, 0]), :] = INF
choices = []
for j in range(m):
    A = P[j]; B = P[(j+1) % m]
    Snew = np.full((G, G), INF); arg = np.zeros((G, G), dtype=np.int8)
    for off in range(-W, W+1):
        # transition k -> k' = k + off
        c = np.full(G, INF)
        if off >= 0:
            diff = B[off:] - A[:G-off]; cc = np.sum(diff**2, 1)/(2*h); c[:G-off] = np.where(np.isnan(cc), INF, cc)
            cand = np.full((G, G), INF); cand[:, off:] = S[:, :G-off] + c[None, :G-off]
        else:
            o = -off
            diff = B[:G-o] - A[o:]; cc = np.sum(diff**2, 1)/(2*h); c[o:] = np.where(np.isnan(cc), INF, cc)
            cand = np.full((G, G), INF); cand[:, :G-o] = S[:, o:] + c[None, o:]
        better = cand < Snew
        Snew[better] = cand[better]; arg[better] = off
    choices.append(arg); S = Snew
prof = np.diag(S).copy()
k0 = int(np.argmin(prof))
# backtrack
ks = [k0]; k = k0
for j in range(m-1, -1, -1):
    off = int(choices[j][k0, k]); k = k - off; ks.append(k)
ks = ks[::-1]
assert ks[0] == k0
t = (np.array(ks[:-1]) + lo)*d
Qopt = np.array([P[j, ks[j]] for j in range(m)])
E = np.sum((np.roll(Qopt, -1, 0)-Qopt)**2)/(2*h)
offs = np.array([int(ks[j+1])-int(ks[j]) for j in range(m)])
print(f"m={m} d={d} margin>={margin}: global grid min E = {prof[k0]:.4f} (recomputed {E:.4f}) at t0 = {(k0+lo)*d:.2f}; "
      f"t range [{t.min():.2f},{t.max():.2f}]; max |offset| used {np.abs(offs).max()} of W={W}")
# local minima of the 1-D landscape E*(t0)
fin = np.isfinite(prof)
tt = (np.arange(G)+lo)*d
loc = [i for i in range(1, G-1) if fin[i] and prof[i] <= prof[i-1] and prof[i] <= prof[i+1] and fin[i-1] and fin[i+1]]
print("local minima of E*(t0):", [(round(tt[i], 2), round(prof[i], 4)) for i in loc])
np.savez(os.path.join(here, f"dp_m{m}_d{d}_mg{margin}.npz"), prof=prof, tt=tt, Qopt=Qopt, t=t)
