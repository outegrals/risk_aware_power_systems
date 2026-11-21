"""Support-wide dispatch checks; replay states from the synthetic microgrid study."""
import itertools, json
from pathlib import Path
import numpy as np
import validate_dynamic_microgrids as d
import validate_risk_and_caps as r

out = Path(__file__).resolve().parents[1]/'output/dispatch'
out.mkdir(exist_ok=True)
C = np.array(r.base_covariance())
L = np.linalg.cholesky(C)
z = np.array(list(itertools.product((-1., 1.), repeat=4)))
records=[]
for s in (0., .5, 1., 2.):
    y = np.array([120.,100.,80.,400.])+s*z@L.T
    pg=y[:,3]-y[:,:3].sum(axis=1)
    assert y.min()>0 and pg.min()>0 and pg.max()<500
    cov=(y-y.mean(axis=0)).T@(y-y.mean(axis=0))/16
    assert np.max(abs(cov-s*s*C))<1e-9
    for q in (np.array([50.,50.,50.]),np.array([100.,50.,50.]),np.array([75.,75.,75.])):
        rt=1+.2*pg; da=1+.2*(400-q.sum())
        profit=(da*q+rt[:,None]*(y[:,:3]-q)).mean(axis=0)
        expected=np.array(r.gaussian_revenue_moments([120,100,80],q.tolist(),(s*s*C).tolist())['individual_expected_revenue'])
        assert np.max(abs(profit-expected))<1e-9
        spread_moment=float(np.mean((rt-da)**2))
        expected_R=r.gaussian_revenue_moments([120,100,80],q.tolist(),(s*s*C).tolist())['price_spread_second_moment_R']
        assert abs(spread_moment-expected_R)<1e-9
    records.append(dict(scale=s,min_generation=y[:,:3].min(),min_load=y[:,3].min(),min_pg=pg.min(),max_pg=pg.max()))
# Same asset seeds and physical-state draws as the 360 published trajectories.
mins=np.full(4,np.inf); count=0
for n in d.SIZES:
 for s in d.SCALES:
  for rep in range(d.REPS):
   a=d.assets(n,d.SEED+rep); energy=.5*a['capacity']; rng=np.random.default_rng(d.SEED+10000+rep)
   capacity=5*a['size'].sum()
   for step in range(d.STEPS):
    p=d.conditional_physics(a,energy,step,s)
    mu,caps,_,load=d.inputs_and_constant(p,5/n)
    signed=p['export']-p['import']
    lo=(p['external_cond_mean']-.02*s*p['external_base']-signed.max(axis=1).sum(axis=1)).min()
    hi=(p['external_cond_mean']+.02*s*p['external_base']-signed.min(axis=1).sum(axis=1)).max()
    # Entire bid box, not only the three selected equilibria.
    margins=np.array([load-caps.sum(),capacity-load,lo,capacity-hi])
    assert margins.min()>0
    mins=np.minimum(mins,margins);count+=1
    common=rng.choice(3,p=d.W);local=rng.choice(3,size=n,p=d.W);rng.choice(3,p=d.W)
    energy=p['next_energy'][common,local,np.arange(n)]
report={'bounded_baseline':records,'baseline_capacity_MW':500,'baseline_day_ahead_range_MW':[100,400],
 'direct_deviation_profit_gain_dollars_per_hour':.2*(100*100-50*150),
 'direct_deviation_J_gain_dollars_per_hour':.2*(100*100-50*150)+(5/3)*.04*(100**2-150**2),
 'dynamic_capacity_rule':'5*sum(asset size) MW','intervals_checked':count,
 'dynamic_min_margins_MW':dict(zip(('DA_lower','DA_upper','RT_lower','RT_upper'),mins.tolist())),
 'scope':'Full conditional support at each replayed state and whole day-ahead bid box; not all possible multistep state histories. Gaussian variance results are not transferred to the bounded distribution.'}
(out/'feasibility_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
