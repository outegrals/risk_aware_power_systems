import json
from pathlib import Path
import numpy as np
rng=np.random.default_rng(182026)
errs=[]
for trial in range(200):
 n=int(rng.integers(2,101));a=float(rng.uniform(.01,1));X=float(rng.uniform(10,1000));u=rng.uniform(.1,2,n)*X/n
 targettotal=min(X/2,u.sum());lo=0.;hi=float(u.max())
 for _ in range(100):
  t=(lo+hi)/2
  if np.minimum(u,t).sum()<targettotal:lo=t
  else:hi=t
 target=np.minimum(u,hi);targettotal=target.sum()
 c=2*a*(X-targettotal-target)
 q=rng.uniform(0,1,n)*u;delta=q-target;T=-c*delta+(np.dot(c,delta)-c*delta)/(n-1)
 assert abs(T.sum())<1e-7
 grad=2*a*(X-target.sum()-target)-c
 assert np.max(abs(grad))<1e-9
 i=int(rng.integers(n));dev=float(rng.uniform(0,u[i]));d=dev-target[i]
 gain=2*a*(dev*(X-targettotal-d)-target[i]*(X-targettotal))-c[i]*d
 assert abs(gain+2*a*d*d)<1e-7
 # Ordinary RN scalar root.
 lo=0.;hi=X
 for _ in range(100):
  r=(lo+hi)/2
  if r+np.minimum(u,r).sum()<X:lo=r
  else:hi=r
 rn=np.minimum(u,(lo+hi)/2)
 p=2*a*target*(X-targettotal);p0=2*a*rn*(X-rn.sum());S=float((p-p0).sum())
 assert S>-1e-7
 b=p0-p+S/n
 assert abs(b.sum())<1e-7
 assert np.max(abs(p+b-p0-S/n))<1e-7
 errs.append(abs(gain+2*a*d*d))
report={'seed':182026,'random_capped_cases':200,'max_unilateral_identity_error':max(errs),'checks':'budget balance, target first-order conditions, unilateral payoff identity, RN fallback surplus, fixed participation-transfer balance','proof_scope':'uniqueness follows from analytic negative-definite potential Hessian, not random sampling'}
Path(__file__).with_name('contracted_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
