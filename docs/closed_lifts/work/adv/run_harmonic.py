import numpy as np, sys
from harmonic import *
ss, Qrest = load_rest(500)
F = Fibres(ss)
Erest = energy(F, Qrest)
print('rest-posture lift: E = %.6f' % Erest)
Q, E, el, w = solve(F, Qrest.copy(), verbose=False)
report('harmonic from rest lift', F, Q, Qrest, E, el, w)
np.save('harm_euclid_500.npy', Q)
# multiple starts: move the rest lift along the self-motion by c rad (Euclidean arc length)
for c in [-4.0, -3.0, -2.0, -1.0, 1.0, 2.0, 3.0, 4.0]:
    Q0 = F.flow(Qrest.copy(), c)
    Qc, Ec, elc, wc = solve(F, Q0)
    report('harmonic from rest lift shifted %+.1f' % c, F, Qc, Qrest, Ec, elc, wc)
    np.save('harm_euclid_500_start%+.0f.npy' % c, Qc)
