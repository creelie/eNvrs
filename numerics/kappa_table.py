# Non-rigorous: growth factor kappa(Lambda) of the renormalized inviscid front (fixed point of the step map)
import numpy as np, sys, time
from renorm import fixed
out=open('kappa_table.txt','w')
out.write('# Lambda kappa s_star alpha_star gamma KB KA fixed_diff\n')
for Lam in [float(v) for v in sys.argv[1].split(',')]:
    r=Lam**(2/3)
    KB=int(min(140, max(40, 25*np.log(10)/np.log(r)/2)))
    KA=8
    t0=time.time()
    try:
        w,G,s,it,d=fixed(Lam,KB,KA,0.5,iters=600)
    except Exception as e:
        print(Lam,'fail',e,flush=True); continue
    a=1/3+np.log(G)/(2*np.log(Lam)); g=1/3-np.log(G)/np.log(Lam)
    line=f"{Lam:.4f} {G:.12f} {s:.8f} {a:.6f} {g:.6f} {KB} {KA} {d:.1e}"
    print(line, f"[{time.time()-t0:.0f}s it={it}]", flush=True); out.write(line+'\n'); out.flush()
