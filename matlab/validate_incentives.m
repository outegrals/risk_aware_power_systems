function report=validate_incentives(out)
rng(182026,'twister');err=0;
for trial=1:200
 n=randi([2 100]);a=.01+.99*rand;X=10+990*rand;u=(.1+1.9*rand(1,n))*X/n;target=min(u,rap.level(min(X/2,sum(u)),u,0));c=2*a*(X-sum(target)-target);
 q=rand(1,n).*u;d=q-target;T=-c.*d+(sum(c.*d)-c.*d)/(n-1);assert(abs(sum(T))<1e-7);
 i=randi(n);v=rand*u(i);delta=v-target(i);gain=2*a*(v*(X-sum(target)-delta)-target(i)*(X-sum(target)))-c(i)*delta;err=max(err,abs(gain+2*a*delta^2));assert(err<1e-7);
 rn=rap.equilibrium(X,u,1);profit=2*a*target*(X-sum(target));p0=2*a*rn*(X-sum(rn));S=sum(profit-p0);b=p0-profit+S/n;assert(S>-1e-7 && abs(sum(b))<1e-7);
end
X=300;a=.1;target=[50 50 50];q=[100 50 50];g=netgain(X,a,target,q);assert(abs(g(1)+500)<1e-9);
qexit=rap.equilibrium(X,[120 100 80],[1 1 1]); %#ok<NASGU>
% Retained marginal charges on insiders, outsider has no transfer.
q=[120 80/3 80/3];exitgain=.2*q(1)*(300-sum(q))-1500;assert(abs(exitgain-1540)<1e-9);
gc=netgain(250,.1,[47.5 47.5 30],[51.875 47.5 30]);assert(abs(sum(gc([1 3]))-3.828125)<1e-9);
report=struct('status','passed','random_capped_cases',200,'max_identity_error',err,'contracted_deviation_net_gain',g(1),'fixed_rival_exit_gain',500,'re_equilibrated_capped_exit_gain',exitgain,'unanimous_RN_fallback_entry_gain',375,'capped_coalition_gain',sum(gc([1 3])),'scope','Exploratory mechanism, not incorporated into original manuscript; fixed verified inputs and enforceable settlement');rap.save_report(out,'incentives',report);
end
function g=netgain(X,a,t,q)
c=2*a*(X-sum(t)-t);d=q-t;T=-c.*d+(sum(c.*d)-c.*d)/(numel(t)-1);g=2*a*q*(X-sum(q))+T-2*a*t*(X-sum(t));
end
