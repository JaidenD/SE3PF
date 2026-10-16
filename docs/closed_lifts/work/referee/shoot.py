"""Referee: one lap of the acceleration-level pseudoinverse ODE q'' = J^+(zeta' - Jdot q') (Theorem B(e), K = I),
with horizontal constraint stabilisation (does not act along ker J), from the stored discrete harmonic lifts' data
(q(0), q'(0) by FFT differentiation). A harmonic lift must close up to (Floquet multiplier) x (discretisation error)."""
import numpy as np, rk
from numpy.linalg import norm, pinv
L = rk.L
def dzeta(s, e=1e-5): return (rk.zeta(s+e)-rk.zeta(s-e))/(2*e)
def acc(s, q, v, a=40.0):
    J = rk.jb(q); e = 1e-6
    Jd = (rk.jb(q+e*v)-rk.jb(q-e*v))/(2*e)
    err = -rk.pose_err(q, s)
    return pinv(J)@(dzeta(s) - Jd@v + 2*a*(rk.zeta(s)-J@v) + a*a*err)
def lap(q, v, n=4000):
    h = L/n; y = np.concatenate([q, v])
    f = lambda s, y: np.concatenate([y[7:], acc(s, y[:7], y[7:])])
    for k in range(n):
        s = k*h; k1 = f(s, y); k2 = f(s+h/2, y+h/2*k1); k3 = f(s+h/2, y+h/2*k2); k4 = f(s+h, y+h*k3); y = y+h/6*(k1+2*k2+2*k3+k4)
    return y[:7], y[7:]
def fftd(Q):
    N = len(Q); w = 2*np.pi*np.fft.fftfreq(N, d=L/N); w[N//2] = 0
    return np.real(np.fft.ifft(1j*w[:, None]*np.fft.fft(Q, axis=0), axis=0))
for name, f in [("A (1000)", "../harmonic_full_1000.npy"), ("B (1000)", None), ("rest lift (1000)", "../rest_lift_1000.npy")]:
    Q = np.load(f) if f else np.load("../num/data/far_I_start-2.npz")["q"]
    V = fftd(Q); qL, vL = lap(Q[0], V[0])
    print(f"{name:18s}: |q(L)-q(0)| = {norm(qL-Q[0]):.2e}, |q'(L)-q'(0)| = {norm(vL-V[0]):.2e}", flush=True)
