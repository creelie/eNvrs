import numpy as np, sys, time
from scipy.integrate import solve_ivp
from scipy.sparse import diags
# critical variables x_n = N_n^{1-2a} u_n, N_n = Lam^n:
#   x_n' = N_n^{2a} (A x_{n-1}^2 - nu x_n - B x_n x_{n+1}),  A = Lam^{1-4a}, B = Lam^{2a-1}, x_{-1} = 0
def xmodel(Lam, a, M, nu=1.0):
    n=np.arange(M+1); R=Lam**(2*a*n); A=Lam**(1-4*a); B=Lam**(2*a-1)
    def f(t,x):
        xm=np.concatenate(([0.0],x[:-1])); xp=np.concatenate((x[1:],[0.0]))
        return R*(A*xm**2-nu*x-B*x*xp)
    def jac(t,x):
        xp=np.concatenate((x[1:],[0.0]))
        return diags([2*R[1:]*A*x[:-1], -R*(nu+B*xp), -R[:-1]*B*x[:-1]],[-1,0,1],format='csc')
    return f,jac
def evolve(Lam, a, M, x0, T, nu=1.0, rtol=1e-9, atol=1e-12, stop_top=1e-6, maxr=5000):
    f,jac=xmodel(Lam,a,M,nu); x=x0.copy(); t=0.0; supx=x.copy(); k=0; first=np.full(M+1,np.nan)
    while t<T and k<maxr:
        sol=solve_ivp(f,(0,T-t),x,method='BDF',jac=jac,rtol=rtol,atol=atol)
        X=sol.y; supx=np.maximum(supx,X.max(1))
        for j in range(len(sol.t)):
            hit=(X[:,j]>1.0)&np.isnan(first); first[hit]=t+sol.t[j]
        x=np.maximum(X[:,-1],0); t+=sol.t[-1]; k+=1
        if sol.status==0 or x[-1]>stop_top: break
    return t,x,supx,first,k
if __name__=='__main__':
    amp=float(sys.argv[1]); Lam=float(sys.argv[2]); alphas=[float(v) for v in sys.argv[3].split(',')]; M=int(sys.argv[4]); T=float(sys.argv[5]) if len(sys.argv)>5 else 5.0
    for a in alphas:
        t0=time.time(); x0=np.zeros(M+1); x0[0]=amp
        t,x,s,first,k=evolve(Lam,a,M,x0,T)
        act=np.nonzero(s>1e-6)[0].max()
        q=[int(act*f) for f in (0.25,0.5,0.75,0.95)]
        print(f"Lam={Lam} a={a:.4f} amp={amp} M={M}: t_end={t:.6g} restarts={k} active to {act}; sup x at {q}: {[f'{s[i]:.4g}' for i in q]}  [{time.time()-t0:.0f}s]", flush=True)
