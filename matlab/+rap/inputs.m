function [mu,caps,constant,load]=inputs(p,a2)
mu=rap.expectation(p.export);n=numel(mu);caps=zeros(1,n);w=[.2 .6 .2];weights=reshape(w'*w,9,1);
for i=1:n
 [v,ix]=sort(reshape(p.export(:,:,i),9,1));j=find(cumsum(weights(ix))>=.1,1);caps(i)=min(mu(i),v(j));
end
signed=p.export-p.import;cm=sum(signed.*reshape(w,1,3,1),2);cv=sum((signed-cm).^2.*reshape(w,1,3,1),2);
r=p.external_cond_mean-reshape(sum(cm,3),3,1);rt=20*rap.expectation(signed)+2*a2*reshape(sum((r.*cm-cv).*reshape(w,3,1,1),1),1,[]);
constant=rt-5*rap.expectation(p.charge+p.discharge);load=p.external_base+sum(rap.expectation(p.import));
assert(max(rap.expectation(p.export<reshape(caps,1,1,n)-1e-12))<=.1+1e-12);
end
