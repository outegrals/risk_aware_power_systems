function report=validate_dispatch_feasibility(out)
C=rap.covariance();L=chol(C,'lower');z=2*rap.product_states(2,4)-1;records=struct([]);i=0;
for s=[0 .5 1 2]
 y=[120 100 80 400]+s*z*L';pg=y(:,4)-sum(y(:,1:3),2);assert(min(y,[],'all')>0 && min(pg)>0 && max(pg)<500);assert(max(abs(cov(y,1)-s*s*C),[],'all')<1e-9);
 for q=[50 100 75;50 50 75;50 50 75]
  q=q';rt=1+.2*pg;da=1+.2*(400-sum(q));v=rap.moments([120 100 80],q,s*s*C);assert(max(abs(mean(da*q+rt.*(y(:,1:3)-q))-v.mean))<1e-9);assert(abs(mean((rt-da).^2)-v.R)<1e-9);
 end
 i=i+1;records=rap.append(records,struct('scale',s,'min_generation',min(y(:,1:3),[],'all'),'min_load',min(y(:,4)),'min_pg',min(pg),'max_pg',max(pg)));
end
f=fullfile(out,'dynamic_microgrids.json');assert(isfile(f),'Run dynamic microgrids first');d=jsondecode(fileread(f));m=min(reshape([d.runs.dispatch_margins_MW],4,[])',[],1); assert(numel(m)==4 && all(m>0));
report=struct('status','passed','bounded_baseline',records,'dynamic_margins_MW',m,'intervals_checked',d.intervals,'scope','Whole bid box and conditional support at replayed states, not every possible history');rap.save_report(out,'dispatch_feasibility',report);
end
