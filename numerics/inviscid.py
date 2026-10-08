import numpy as np, sys, time
from scipy.integrate import solve_ivp
from scipy.sparse import diags
# inviscid critical lattice: y_n' = r^n (y_{n-1}^2 - y_n y_{n+1}), r = Lam^{2/3} (rate ratio)
def front_growth(Lam, M=400, T=1e3):
    r=Lam**(2/3); R=r**np.arange(M+1)
    def f(t,y):
        ym=np.concatenate(([0.0],y[:-1])); yp=np.concatenate((y[1:],[0.0]))
        return R*(ym**2-y*yp)
    def jac(t,y):
        yp=np.concatenate((y[1:],[0.0]))
        return diags([2*R[1:]*y[:-1], -R*yp, -R[:-1]*y[:-1]],[-1,0,1],format='csc')
    y=np.zeros(M+1); y[0]=1.0; t=0.0; s=y.copy(); k=0
    while t<T and k<5000:
        sol=solve_ivp(f,(0,T-t),y,method='BDF',jac=jac,rtol=1e-10,atol=1e-14)
        s=np.maximum(s,sol.y.max(1)); y=np.maximum(sol.y[:,-1],0); t+=sol.t[-1]; k+=1
        if sol.status==0 or y[-1]>1e-8: break
    act=np.nonzero(s>1e-8)[0].max()
    i1,i2=int(0.4*act),int(0.8*act)
    g=np.exp((np.log(s[i2])-np.log(s[i1]))/(i2-i1))
    return g, act, s
if __name__=='__main__':
    for Lam in [float(v) for v in sys.argv[1].split(',')]:
        t0=time.time(); g,act,s=front_growth(Lam, M=int(sys.argv[2]))
        print(f"Lam={Lam}: front growth per shell kappa={g:.6f} (active {act}); predicted alpha*=1/3+log(kappa)/(2 log Lam)={1/3+np.log(g)/(2*np.log(Lam)):.4f}; sup y at 10,50,100,200: {[f'{s[i]:.4g}' for i in (10,50,100,200) if i<len(s)]} [{time.time()-t0:.0f}s]", flush=True)
