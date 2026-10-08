import numpy as np, sys, time
from inviscid import front_growth
for Lam in [float(v) for v in sys.argv[1].split(',')]:
    r=Lam**(2/3); M=int(min(400, 60*np.log(10)/np.log(r)))
    t0=time.time(); g,act,s=front_growth(Lam, M=M)
    print(f"Lam={Lam}: M={M} kappa={g:.6f} act={act} alpha*={1/3+np.log(g)/(2*np.log(Lam)):.4f} gamma={1/3-np.log(g)/np.log(Lam):.4f} [{time.time()-t0:.0f}s]", flush=True)
