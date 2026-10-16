# Is the harmonic lift (K = I) a periodic orbit of the acceleration-level pseudoinverse
#   q'' = J^+ (zeta_gamma' - Jdot q' + kd (zeta_gamma - J q') + kp e) ?
# (the stabilizing terms lie in range J^+ = horizontal, so the vertical part of q'' is still zero)
import numpy as np
from pandageo import *
from lifts import *
path = ConingCircle(); L = path.L
Qh = np.load('harmonic_full_1000.npy'); m = Qh.shape[0]; h = L/m
q0 = Qh[0].copy()
dq0 = (-Qh[2] + 8*Qh[1] - 8*Qh[-1] + Qh[-2])/(12*h)
def acc(q, dq, s, kp=400.0, kd=40.0):
    e, J = pose_error(q, path, s)
    eps = 1e-6
    Jd = (fk(q+eps*dq)[2] - fk(q-eps*dq)[2])/(2*eps)
    rhs = path.acceleration(s) - Jd@dq + kd*(path.velocity(s) - J@dq) + kp*e
    return pinv(J)@rhs
def shoot(q, dq, n=4000):
    hh = L/n; s = 0.0
    for i in range(n):
        f = lambda y, s: np.concatenate([y[7:], acc(y[:7], y[7:], s)])
        y = np.concatenate([q, dq])
        k1=f(y,s); k2=f(y+0.5*hh*k1,s+0.5*hh); k3=f(y+0.5*hh*k2,s+0.5*hh); k4=f(y+hh*k3,s+hh)
        y = y + (hh/6)*(k1+2*k2+2*k3+k4); s += hh
        q, dq = y[:7], y[7:]
        if i in (n//4, n//2, 3*n//4):
            j = int(round(s/h)) % m
            print(f'   s = {s:.4f}: |q - q_harm(s)| = {np.linalg.norm(q - Qh[j]):.2e}')
    return q, dq
qL, dqL = shoot(q0, dq0)
print('one lap of the acceleration-level pseudoinverse from (q*(0), q*\'(0)):')
print('  |q(L) - q*(0)| =', np.linalg.norm(qL-q0), '  |q\'(L) - q*\'(0)| =', np.linalg.norm(dqL-dq0))
# compare: same from the rest-posture lift
Qr = np.load('rest_1000.npy'); dqr = (-Qr[2] + 8*Qr[1] - 8*Qr[-1] + Qr[-2])/(12*h)
qL, dqL = shoot(Qr[0].copy(), dqr)
print('from the rest-posture lift data: |q(L) - q_rest(0)| =', np.linalg.norm(qL-Qr[0]), ' |q\'(L)-q\'(0)| =', np.linalg.norm(dqL-dqr))
