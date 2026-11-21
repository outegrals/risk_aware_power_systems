function r=network(g,k,loss,onepeer,episodes,slot)
n=g.m-1;
if onepeer,batches=g.colors;if mod(n,2)==0,setup=[n/2 n/2];else,setup=[floor(n/2) floor(n/2) 1];end
else,batches=2*size(g.edges,1);setup=n;end
sizes=unique([batches setup]);slots=zeros(episodes,1);expected=0;tail=0;
for messages=sizes
 count=k*sum(batches==messages)+sum(setup==messages);
 if loss==0,slots=slots+count;expected=expected+count;continue;end
 amax=max(1,ceil(log(1e-13/(messages*count))/log(loss)));cdf=exp(messages*log1p(-loss.^(1:amax)));pmf=diff([0 cdf]);pmf(end)=pmf(end)+1-sum(pmf);
 remain=repmat(count,episodes,1);mass=1;
 for j=1:numel(pmf)-1
  draws=rap.binomial(remain,min(1,max(0,pmf(j)/mass)));slots=slots+j*draws;remain=remain-draws;mass=mass-pmf(j);
 end
 slots=slots+numel(pmf)*remain;expected=expected+count*(1+sum(-expm1(messages*log1p(-loss.^(1:amax+9)))));tail=tail+messages*loss^amax*count;
end
seconds=slots*slot/1000;deadlines=[1 10 60 300 900];f=mean(seconds<=deadlines,1);z=1.959963984540054;center=(f+z*z/(2*episodes))/(1+z*z/episodes);radius=z*sqrt(f.*(1-f)/episodes+z*z/(4*episodes^2))/(1+z*z/episodes);
r=struct('one_peer',onepeer,'loss',loss,'slot_ms',slot,'episodes',episodes,'expected_seconds',expected*slot/1000,'mean_seconds',mean(seconds),'p50_seconds',median(seconds),'p95_seconds',quantile(seconds,.95),'deadlines',deadlines,'met_fraction',f,'wilson_lower',max(0,center-radius),'wilson_upper',min(1,center+radius),'tail_bound',tail,'expected_attempts',(2*size(g.edges,1)*k+n)/(1-loss));
end
