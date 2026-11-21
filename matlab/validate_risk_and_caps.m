function report=validate_risk_and_caps(out)
% Native moment quadrature, constrained best responses and coefficient retuning.
x=[120 100 80];C=rap.covariance();L=chol(C,'lower');states=rap.product_states(3,4)+1;
z=[-sqrt(3) 0 sqrt(3)];w=[1/6 2/3 1/6];weights=prod(w(states),2);records=struct([]);r=0;
for scale=[0 .5 1 2]
 y=[x 400]+scale*z(states)*L';rt=1+.2*(y(:,4)-sum(y(:,1:3),2));
 for q=[75 50;75 50;75 50]
  q=q';v=rap.moments(x,q,scale^2*C);da=1+.2*(400-sum(q));rev=da*q+rt.*(y(:,1:3)-q);m=weights'*rev;varr=weights'*(rev-m).^2;R=weights'*(rt-da).^2;
  err=[max(abs(m-v.mean)),max(abs(varr-v.variance)),abs(R-v.R)];assert(all(err<[1e-8 1e-7 1e-9]));
  k=q(1)==75;if k,k=1;else,k=1/3;end
  for i=1:3
   vals=zeros(1,3);
   for j=1:3
    trial=q;trial(i)=(j-1)*x(i)/2;price=1+.2*(400-sum(trial));rr=price*trial+rt.*(y(:,1:3)-trial);
    vals(j)=weights'*rr(:,i)-(k-1)/.4*(weights'*(rt-price).^2);
   end
   h=x(i)/2;A=(vals(3)-2*vals(2)+vals(1))/(2*h*h);B=(vals(2)-vals(1)-A*h*h)/h;
   assert(abs(min(x(i),max(0,-B/(2*A)))-q(i))<1e-8);
  end
  r=r+1;records=rap.append(records,struct('scale',scale,'q',q,'moments',v,'quadrature_errors',err));
 end
end
xs={[120 100 80],[120 100 30],[210 45 45],[1 18 18 18 45],[120 100 80],[120 100 80]};us=xs;us{5}=[30 20 10];us{6}=[0 0 0];known=struct([]);
for i=1:6
 [q,k,res]=rap.strategies(xs{i},us{i});known=rap.append(known,struct('x',xs{i},'caps',us{i},'q',q,'k',k,'residual',res));
end
assert(max(abs(known(2).q(3,:)-[47.5 47.5 30]))<1e-10);
assert(max(abs(known(4).q(3,:)-[1 12.25 12.25 12.25 12.25]))<1e-10);
rng(20261117,'twister');maxerr=0;
for n=[2 3 10 50 100]
 for trial=1:20
  x=exp(log(1)+rand(1,n)*log(200));u=x.*rand(1,n)*1.5;[q,k,res]=rap.strategies(x,u);maxerr=max(maxerr,res);
  for j=1:3,assert(max(abs(q(j,:)-rap.equilibrium(sum(x),u,k(j))))<1e-8);end
 end
end
report=struct('status','passed','uncertainty',records,'known_caps',known,'random_cases',100,'max_br_residual',maxerr,'rng','MATLAB twister; randomized cases differ from Python');
rap.save_report(out,'risk_and_caps',report);
end
