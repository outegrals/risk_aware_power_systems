import json

def payoff(X,a,target,q):
 n=len(q);Q=sum(q);Qt=sum(target)
 c=[2*a*(X-Qt-t) for t in target]
 charges=[c[i]*(q[i]-target[i]) for i in range(n)]
 transfers=[-charges[i]+sum(charges[j] for j in range(n) if j!=i)/(n-1) for i in range(n)]
 return [2*a*q[i]*(X-Q)+transfers[i] for i in range(n)],transfers
out={}
for name,X,target,q in [('symmetric',300,[50,50,50],[62.5,62.5,50]),('capped',250,[47.5,47.5,30],[51.875,47.5,30])]:
 old,_=payoff(X,.1,target,target);new,T=payoff(X,.1,target,q)
 out[name]={'net_gains':[new[i]-old[i] for i in range(3)],'transfers':T,'budget_sum':sum(T)}
assert out['symmetric']['net_gains']==[-62.5,-62.5,0]
assert abs(out['capped']['net_gains'][0]+out['capped']['net_gains'][2]-3.828125)<1e-9
print(json.dumps(out,indent=2))
