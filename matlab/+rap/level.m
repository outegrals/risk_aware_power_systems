function t=level(total,u,inverse_k)
u=sort(u(:)'); n=numel(u); prefix=[0 cumsum(u)];
if inverse_k==0 && total>=sum(u)-1e-12, t=max(u); return; end
for j=0:n
 den=inverse_k+n-j;
 if den<=0, continue; end
 t=max(0,(total-prefix(j+1))/den);
 if j==n || t<=u(j+1)+1e-12, return; end
end
error('No root');
end
