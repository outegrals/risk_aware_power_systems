function z=step(z,g)
n=g.m-1;
switch g.kind
 case 'star',z=[(n*z(1:n,:)+z(end,:))/(n+1);sum(z,1)/(n+1)];
 case 'complete',z=repmat(mean(z,1),g.m,1);
 otherwise
  v=z;for s=g.shifts,v=v+circshift(z,s,1);end;z=v/g.denominator;
end
end
