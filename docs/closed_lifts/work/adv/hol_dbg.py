import numpy as np
exec(open('holonomy.py').read().split("# locate")[0].replace("M = 140", "M = 140"))
from selfmotion import wrap
dd = []
for q in Q:
    dd.append(np.linalg.norm(wrap(QQ - q), axis=1).min())
dd = np.array(dd); j = int(np.argmax(dd)); print('worst', j, dd[j], 'theta0', th0[j])
print('q0', QQ[idx[j]]); print('q1', Q[j]); print('wrapped diff to start', wrap(Q[j]-QQ[idx[j]]))
