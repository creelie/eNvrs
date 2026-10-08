# Builds the centre and the box radii of one case, then runs the verifier on it.
# usage (from cap/): python3 make_case.py Lnum Lden KB KA aHiNum aHiDen delta theta_contr
#   centre  : fixed point of the inviscid step map (numerics/renorm.py), not rigorous;
#   radii   : delta * (I - |J|/theta_contr)^{-1} b, J the step Jacobian at the centre (from the
#             verifier's point mode), b = 1 behind the front and max(0.1 w_k, AFLOOR) ahead of it.
# The verifier then proves the step for the whole box; nothing computed here needs to be rigorous.
# The radii files of Lambda = 13/10 come from an earlier run with b_k = max(10^ceil(log10 w_k), 1e-11)
# ahead of the front; the other four cases were made by this script (AFLOOR=1e-11 for Lambda = 9/5).
import sys, subprocess, numpy as np
sys.path.insert(0,'../numerics')
from renorm import fixed
Ln,Ld,KB,KA=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]); an,ad=sys.argv[5],sys.argv[6]
delta=float(sys.argv[7]); th=float(sys.argv[8]); Lam=Ln/Ld
tag=f"L{Ln}_{Ld}_K{KB}_{KA}"
w,G,s,it,d=fixed(Lam,KB,KA,0.5,iters=400)
ks=np.arange(-KB,KA+1)
np.savetxt(f'center_{tag}.txt', np.c_[ks,w], fmt=['%d','%.17g'])
print(tag,'nonrigorous kappa',G,'s*',s,'fixed diff',d, flush=True)
subprocess.run(['./verify_step',str(KB),str(KA),str(Ln),str(Ld),'0','1',an,ad,'1e-12',f'center_{tag}.txt','none',f'pt_{tag}','1.0','1e-31'],check=True)
J=np.loadtxt(f'pt_{tag}_Jcenter.txt')
keep=(ks!=0)&(ks!=1); A=np.abs(J)[np.ix_(keep,keep)]; kk=list(ks[keep])
print('rho(|J|)=',max(abs(np.linalg.eigvals(A))))
wa={k:abs(w[KB+k]) for k in range(2,KA+1)}
import os
fl=float(os.environ.get('AFLOOR','1e-13'))
b=np.array([1.0 if k<0 else max(wa[k]*1e-1, fl) if k<KA else fl for k in kk])
# ahead weights: keep at least ~ value*0.1 for k>=2 but never below 1e-13
rho=np.linalg.solve(np.eye(len(kk))-A/th,b)
R=delta*rho
R[kk.index(KA)]=1e-20
with open(f'radii_{tag}.txt','w') as f:
    for k,v in zip(kk,R): f.write(f"{k} {v:.17g}\n")
out=subprocess.run(['./verify_step',str(KB),str(KA),str(Ln),str(Ld),'0','1',an,ad,'1e-12',f'center_{tag}.txt',f'radii_{tag}.txt',f'box_{tag}','1.0','1e-31'],capture_output=True,text=True)
print(out.stdout.strip())
print(open(f'box_{tag}.log').read())
