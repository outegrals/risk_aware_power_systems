function [q,k,res]=strategies(mu,u)
X=sum(mu);n=numel(mu);target=min(X/2,sum(u));t=rap.level(target,u,0);
k=[1;1/n;t/(X-target)];
q=[min(u,rap.level(X,u,1));min(u,rap.level(X,u,n));min(u,t)];
br=min(u,max(0,k./(1+k).*(X-sum(q,2)+q)));
res=max(abs(br-q),[],'all')/max(1,X);
assert(res<1e-10 && abs(sum(q(3,:))-target)<1e-10*max(1,X));
end
