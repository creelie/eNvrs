import numpy as np, sys
from scipy.integrate import solve_ivp
# normalized front system, window k=-KB..KA (w_0=1 fixed), section w_1=c
def make(Lam, KB, KA, eps0=0.0, alpha=1/3, wdeep=0.0):
    r=Lam**(2/3); ks=np.arange(-KB,KA+1); rk=r**ks; i0=KB
    vk=Lam**(2*alpha*ks)-1.0
    def f(s,z):
        w=z[:-2]; eps=z[-2]  # z[-1] = log growth
        wm=np.concatenate(([wdeep],w[:-1])); wp=np.concatenate((w[1:],[0.0]))
        g=w[i0-1]**2-w[i0+1]
        dw=rk*(wm**2-w*wp)-w*g-eps*vk*w
        dw[i0]=0.0
        return np.concatenate((dw,[-eps*(g-eps), g-eps]))
    return f, ks, i0
def step(Lam,KB,KA,w,c,eps=0.0,alpha=1/3,rtol=1e-13,atol=1e-16):
    f,ks,i0=make(Lam,KB,KA,eps,alpha)
    z0=np.concatenate((w,[eps,0.0]))
    def ev(s,z): return z[i0+2]-c*z[i0+1]
    ev.terminal=True; ev.direction=1
    sol=solve_ivp(f,(0,100),z0,method='DOP853',rtol=rtol,atol=atol,events=ev)
    z=sol.y_events[0][0]; s=sol.t_events[0][0]
    w1=z[i0+1]; wn=np.concatenate((z[1:len(w)]/w1,[0.0]))
    G=w1*np.exp(z[-1])
    return wn, G, s, z[-2]*1.0
def fixed(Lam,KB,KA,c,iters=200,w=None,tol=1e-14):
    if w is None:
        ks=np.arange(-KB,KA+1); w=np.where(ks<=0, 1.02**(ks*1.0), 0.0); w[KB+1]=c
    for it in range(iters):
        wn,G,s,_=step(Lam,KB,KA,w,c)
        d=np.abs(wn-w).max(); w=wn
        if d<tol: break
    return w,G,s,it,d
if __name__=='__main__':
    Lam=float(sys.argv[1]); c=float(sys.argv[2]); KA=int(sys.argv[3])
    for KB in [int(v) for v in sys.argv[4].split(',')]:
        w,G,s,it,d=fixed(Lam,KB,KA,c)
        print(f"Lam={Lam} c={c} KB={KB} KA={KA}: kappa={G:.12f} step s*={s:.6f} iters={it} last diff={d:.1e}")
        print("  w[-6..KA]:", " ".join(f"{v:.6g}" for v in w[KB-6:]))
        np.save(f"fix_{Lam}_{c}_{KB}_{KA}.npy", w)
