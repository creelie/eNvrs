# Front profiles for the 3D figure: full-lattice front iteration (Lambda = 13/10, alpha = 3/8, eps0 = 1e-12)
# from the centre of the verified box; records log10 y_m at the start of every 10th step.
import numpy as np
from scipy.integrate import solve_ivp
Lam=13/10; r=Lam**(2/3); c=0.5; KB=60; KA=6; AH=12; alpha=0.375; eps=1e-12
cen=np.loadtxt('../cap/cases/center_L13_10_K60_6.txt')
n=KB; y=np.zeros(n+AH+1)
for k,v in cen: y[n+int(k)]=v
logY=0.0; tau=0.0; rows=[]
for step in range(101):
    if step%10==0:
        for m in range(max(0,n-60), n+9):
            v=np.exp(logY)*y[m] if m < len(y) else 0.0
            rows.append((step, m, np.log10(max(v,1e-14)), tau))
    M=len(y); idx=np.arange(M)-n; R=r**idx.astype(float); E=eps*Lam**(2*alpha*idx)
    def f(s,Y):
        Ym=np.concatenate(([0.0],Y[:-1])); Yp=np.concatenate((Y[1:],[0.0]))
        return R*(Ym**2-Y*Yp)-E*Y
    def ev(s,Y): return Y[n+2]-c*Y[n+1]
    ev.terminal=True; ev.direction=1
    sol=solve_ivp(f,(0,5),y,method='DOP853',rtol=1e-12,atol=1e-30,events=ev)
    Y=sol.y_events[0][0]; s=sol.t_events[0][0]; G=Y[n+1]
    tau+=s*np.exp(-(n-KB)*np.log(r)-logY)   # in units of 1/(r^{n0} Y_{n0})
    logY+=np.log(G); y=np.concatenate((Y/G,[0.0])); eps*=Lam**(2*alpha)/(r*G); n+=1
np.savetxt('../paper/figures/data/profiles.dat', np.array(rows), fmt=['%d','%d','%.6f','%.10f'], header='step shell log10y tau')
print('final logY/log(10)=',logY/np.log(10),' tau=',tau)
# per-step files read by paper/figures/fig_front.tex: log10 y clipped at -1.5, cut after the first clipped shell
P = np.array(rows)
for j in range(0, 101, 10):
    with open(f'../paper/figures/data/front_{j}.dat', 'w') as fo:
        fo.write('m n z\n')
        for st, m, z, _ in P[P[:, 0] == j]:
            zc = max(z, -1.5)
            fo.write(f'{int(m)} {int(st)} {zc:.5f}\n')
            if zc == -1.5:
                break
