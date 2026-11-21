function v=moments(x,q,C)
a=.1;X=sum(x);Q=sum(q);V=C(4,4)-2*sum(C(1:3,4))+sum(C(1:3,1:3),'all');vp=4*a*a*V;mp=1+2*a*(400-X);
c=2*a*(C(1:3,4)'-sum(C(1:3,1:3),2)');exposure=x-q;
v.mean=mp*x+c+2*a*q*(X-Q);
v.variance=mp^2*diag(C(1:3,1:3))'+exposure.^2*vp+2*mp*exposure.*c+vp*diag(C(1:3,1:3))'+c.^2;
v.R=vp+4*a*a*(X-Q)^2;v.price_variance=vp;
end
