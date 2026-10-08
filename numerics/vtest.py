import numpy as np, sys, time
from xmodel import evolve
# start from the renormalized front profile (inviscid fixed point) placed at shell n0, amplitude amp
Lam=float(sys.argv[1]); amp=float(sys.argv[2]); M=int(sys.argv[3]); alphas=[float(a) for a in sys.argv[4].split(',')]
w=np.load(f"fix_1.3_0.5_40_8.npy"); KB=40
for a in alphas:
    th=Lam**(2/3-2*a); x0=np.zeros(M+1); n0=KB
    for i,v in enumerate(w):
        x0[i]=amp*th**i*v
    t0=time.time(); t,x,s,first,k=evolve(Lam,a,M,x0,T=50.0)
    act=np.nonzero(s>1e-6*amp)[0].max()
    idx=[n0,n0+50,n0+100,n0+200,n0+400,n0+600,min(act,M-5)]
    print(f"a={a}: t_end={t:.4g} active={act}; sup x/amp at {idx}: {[f'{s[i]/amp:.4g}' for i in idx if i<=M]} [{time.time()-t0:.0f}s]",flush=True)
