import numpy as np
def logm_so(R):
    """principal real logarithm of a rotation matrix via real Schur-like eigen decomposition (no 180-degree angles)."""
    w, V = np.linalg.eig(R)
    Lg = V @ np.diag(np.log(w.astype(complex))) @ np.linalg.inv(V)
    Lg = np.real(Lg); return 0.5*(Lg-Lg.T)
def expm(A, terms=40):
    # scaling and squaring Taylor
    nrm = np.abs(A).sum(); k = max(0, int(np.ceil(np.log2(nrm+1e-300)))+1) if nrm > 0.5 else 0
    B = A/2**k; E = np.eye(len(A)); T = np.eye(len(A))
    for i in range(1, terms): T = T@B/i; E = E+T
    for _ in range(k): E = E@E
    return E
