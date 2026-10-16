import numpy as np
from kin import QR, QMIN, QMAX, Wd
from harmonic import Fibres, load_rest, grad_hess, energy
print('--- (i) Morse count of f = |q - q_r|^2/2 on traced self-motion circles; (ii) joint winding of the circle')
for k in [0, 50, 100, 150, 200, 238, 250, 300, 350, 400, 450]:
    d = np.load('sm_%d.npz' % k); TH, QQ = d['TH'], d['QQ']
    # unwrapped joint path along the circle: QQ is continuous already (no wrapping in trace)
    wind = (QQ[-1] - QQ[0]) / (2 * np.pi)
    # f on the circle with wrapping to the representative nearest q_r (joints on T^7)
    Qw = QR + (QQ - QR + np.pi) % (2 * np.pi) - np.pi
    f = 0.5 * np.sum((Qw - QR) ** 2, 1)
    n = len(f); mins = [i for i in range(n) if f[i] < f[i - 1] and f[i] <= f[(i + 1) % n]]
    maxs = [i for i in range(n) if f[i] > f[i - 1] and f[i] >= f[(i + 1) % n]]
    feas = np.all((Qw > QMIN) & (Qw < QMAX), 1)
    print(f"k={k:3d}: joint winding over circle {np.round(wind, 3)}, local minima of f at theta {[round(TH[i], 3) for i in mins]} "
          f"(f={[round(f[i], 3) for i in mins]}, feasible {[bool(feas[i]) for i in mins]}), #maxima {len(maxs)}")
print('--- (v) Theorem A: connection form A = N^T W N\' along the paper frame (should be constant = -Phi/L)')
fr = np.genfromtxt('../num/data/frame_500.csv', delimiter=',', names=True)
Ns = np.stack([fr['N%d' % i] for i in range(30)], 1).reshape(-1, 5, 6).transpose(0, 2, 1)   # column-major 6x5
dNs = np.stack([fr['dN%d' % i] for i in range(30)], 1).reshape(-1, 5, 6).transpose(0, 2, 1)
W = np.diag(Wd)
A = np.einsum('kia,ij,kjb->kab', Ns, W, dNs)
print('max |A - A^T|_skewness check', np.abs(A + A.transpose(0, 2, 1)).max(), ' variation max_k |A_k - A_0|', np.abs(A - A[0]).max())
ev = np.linalg.eigvals(A[0] * 0.7700876453125667)
print('eigenvalues of L*A (= -Phi): ', np.round(np.sort_complex(ev), 8))
print('--- (vi)/(vii) lifts: max|q\'|, Floquet multipliers of the discrete Jacobi operator (minimizer => hyperbolic)')
ss, Qrest = load_rest(500); F = Fibres(ss)
mp = np.load('midpoint_500.npz')
lifts = {'rest': Qrest, 'midpoint': mp['Qm'], 'harmonic A (E=6.669)': np.load('harm_euclid_500.npy'),
         'harmonic B (E=5.128)': np.load('harm_euclid_500_start-3.npy'), 'harmonic C (E=6.784)': np.load('harm_euclid_500_start+3.npy')}
for name, Q in lifts.items():
    dq = (np.roll(Q, -1, 0) - np.roll(Q, 1, 0)) / (2 * F.h)
    line = f"{name:22s}: E={energy(F, Q):.4f}  max|q'|={np.linalg.norm(dq, axis=1).max():.3f}"
    if name.startswith('harmonic'):
        g, H, X = grad_hess(F, Q)
        N = F.N; d = np.diag(H); o = np.array([H[k, (k + 1) % N] for k in range(N)])
        M = np.eye(2)
        for k in range(N):   # o_{k-1} y_{k-1} + d_k y_k + o_k y_{k+1} = 0
            T = np.array([[-d[k] / o[k], -o[k - 1] / o[k]], [1.0, 0.0]])
            M = T @ M
        mu = np.linalg.eigvals(M)
        line += f"  monodromy trace {np.trace(M):.4e}, multipliers {mu.real.max():.4e}, {mu.real.min():.4e}"
    print(line)
