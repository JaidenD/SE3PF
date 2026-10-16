import numpy as np
from scipy.linalg import logm, expm
from pandageo import *
path = ConingCircle(); L = path.L; W = np.diag(path.Wd)
n = 2000; h = L/n
z0 = path.velocity(0.0)
# Gram-Schmidt in W on zeta(0), e1..e6
basis=[]
for v in [z0]+[np.eye(6)[i] for i in range(6)]:
    w = v.copy()
    for b in basis: w = w - (b@W@w)*b
    nw = np.sqrt(w@W@w)
    if nw > 1e-8: basis.append(w/nw)
N0 = np.array(basis[1:]).T
def f(N, s):
    z = path.velocity(s); dz = path.acceleration(s)
    return -np.outer(z, dz@W@N)
Ns=[N0.copy()]; N=N0.copy(); s=0
for i in range(n):
    k1=f(N,s); k2=f(N+0.5*h*k1,s+0.5*h); k3=f(N+0.5*h*k2,s+0.5*h); k4=f(N+h*k3,s+h)
    N=N+(h/6)*(k1+2*k2+2*k3+k4); s+=h; Ns.append(N.copy())
Ns=np.array(Ns); sg=np.arange(n+1)*h
H = N0.T@W@Ns[-1]
U,S,Vt=np.linalg.svd(H); Hr=U@Vt
Phi = np.real(logm(Hr)); Phi=0.5*(Phi-Phi.T)
ev = np.linalg.eigvals(Hr); th = np.sort(np.abs(np.angle(ev)))[::-1]
print('holonomy angles (pairs):', th)
theta = th[[0,2]]
print('N0^T W N0prime (should be 0):', max(np.abs(Ns[k].T@W@f(Ns[k],sg[k])).max() for k in range(0,n+1,100)))
print('|Phi|_F^2/(2L) =', np.sum(Phi**2)/(2*L), '  sum theta^2 / L =', np.sum(theta**2)/L)
# energy of the paper's frame by quadrature, and of perturbed closed frames
def energy(Qfun, dQfun):
    E=0
    for k in range(n+1):
        Q=Qfun(sg[k]); dQ=dQfun(sg[k]); Nk=Ns[k]@Q
        dN = f(Ns[k],sg[k])@Q + Ns[k]@dQ
        A = Nk.T@W@dN
        wgt = 0.5 if k in (0,n) else 1.0
        E += wgt*h*0.5*np.sum(A**2)
    return E
Ep = energy(lambda s: expm(-(s/L)*Phi), lambda s: -(1/L)*Phi@expm(-(s/L)*Phi))
print('energy of N0 exp(-s/L Phi) by quadrature:', Ep)
rng=np.random.default_rng(1)
for trial in range(3):
    X=rng.normal(size=(5,5)); X=X-X.T
    eps=0.05; kfreq=trial+1
    Qf=lambda s: expm(-(s/L)*Phi)@expm(eps*np.sin(2*np.pi*kfreq*s/L)*X)
    dQf=lambda s: (-(1/L)*Phi@expm(-(s/L)*Phi)@expm(eps*np.sin(2*np.pi*kfreq*s/L)*X)
                   + expm(-(s/L)*Phi)@expm(eps*np.sin(2*np.pi*kfreq*s/L)*X)@(eps*2*np.pi*kfreq/L*np.cos(2*np.pi*kfreq*s/L)*X))
    print('  perturbed closed frame', trial, 'energy', energy(Qf,dQf))
# full vs vertical energy identity |N'|_W^2 = |zeta'|_W^2 + |N^T W N'|^2
k=777; Q=expm(-(sg[k]/L)*Phi); Nk=Ns[k]@Q; dN=f(Ns[k],sg[k])@Q+Ns[k]@(-(1/L)*Phi@Q)
dz=path.acceleration(sg[k])
print('identity check:', np.trace(dN.T@W@dN), dz@W@dz + np.sum((Nk.T@W@dN)**2))
