# Independent (non-rigorous) check of the induction: full lattice (all deep shells, 12 ahead shells),
# viscous term included, start from a corner of the box B, iterate many front steps.
import numpy as np, sys
from scipy.integrate import solve_ivp
Lam=13/10; r=Lam**(2/3); c=0.5; KB=60; KA=6; AH=12
alpha=float(sys.argv[1]); eps=float(sys.argv[2]); nsteps=int(sys.argv[3]); sign=int(sys.argv[4])
cen=np.loadtxt('../cap/cases/center_L13_10_K60_6.txt'); rad=dict((int(k),v) for k,v in np.loadtxt('../cap/cases/radii_L13_10_K60_6.txt'))
# state: shells 0..n0+AH, front index n0=KB
n=KB; y=np.zeros(n+AH+1)
for k,v in cen: y[n+int(k)]=v
rng=np.random.default_rng(sign)
for k in range(-KB,KA+1):
    if k in (0,1): continue
    s=rng.choice([-1,1]) if sign>0 else 1
    y[n+k]+= s*0.99*rad[k] if k<KA else 0.99*rad[k]
y[y<0]=0
worst=0; Gs=[]
for step in range(nsteps):
    m=len(y); idx=np.arange(m)-n; R=r**idx.astype(float); E=eps*Lam**(2*alpha*idx)
    def f(s,Y):
        Ym=np.concatenate(([0.0],Y[:-1])); Yp=np.concatenate((Y[1:],[0.0]))
        return R*(Ym**2-Y*Yp)-E*Y
    def ev(s,Y): return Y[n+2]-c*Y[n+1]
    ev.terminal=True; ev.direction=1
    sol=solve_ivp(f,(0,5),y,method='DOP853',rtol=1e-12,atol=1e-30,events=ev)
    Y=sol.y_events[0][0]; G=Y[n+1]; Gs.append(G)
    y=np.concatenate((Y/G,[0.0])); eps=eps*Lam**(2*alpha)/(r*G); n=n+1
    # compare window with box
    for k in range(-KB,KA+1):
        if k in (0,1): continue
        dev=abs(y[n+k]-cen[KB+k,1])/rad[k]; 
        if dev>worst: worst=dev; wk=k
    if step%25==0 or step==nsteps-1:
        print(f"step {step}: G={G:.12f} eps={eps:.3e} worst box ratio so far {worst:.4f} (k={wk}) deep max {y[:n-KB].max() if n>KB else 0:.4f} ahead max {y[n+KA+1:].max():.2e}",flush=True)
print("min G",min(Gs),"max G",max(Gs))
