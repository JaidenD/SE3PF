import numpy as np
from harmonic import Fibres, load_rest, solve
ss, Qrest = load_rest(500); F = Fibres(ss)
for kp in [10.0, 100.0, 1000.0, 10000.0]:
    Q, E, el, w = solve(F, Qrest.copy(), kp=kp)
    dist = np.linalg.norm(Q - Qrest, axis=1).max()
    print(f"k={kp:8.0f}: max|q_k - q*_rest| = {dist:.4e}, k*dist = {kp*dist:.3f}, Hess min eig {w[0]:.3e}, EL {el:.1e}", flush=True)
