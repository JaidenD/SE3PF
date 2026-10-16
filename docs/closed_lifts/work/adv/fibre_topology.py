"""Does the full self-motion fibre (joints on T^7, limits ignored) change topology along the path?
Trace the circle through the rest-lift node with a small step and record min sigma_min and joint winding."""
import numpy as np, sys
from selfmotion import trace, wrap
from kin import fk
rest = np.genfromtxt('../num/data/rest_500.csv', delimiter=',', names=True)
qs = np.stack([rest['q%d' % i] for i in range(7)], 1); ss = rest['s']
ks = [int(a) for a in sys.argv[1:]] if len(sys.argv) > 1 else list(range(150, 205, 5))
for k in ks:
    TH, QQ, SM = trace(qs[k], ss[k], h=0.002)
    wind = np.round((QQ[-1] - QQ[0]) / (2 * np.pi), 3)
    j = int(np.argmin(SM))
    print(f"k={k:3d} s={ss[k]:.4f}: length {TH[-1]:.4f}, winding {wind}, min sigma_min {SM.min():.5f} at theta {TH[j]:.3f}, q there {np.round(QQ[j], 3)}")
    sys.stdout.flush()
