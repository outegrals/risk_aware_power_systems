function report=audit_dynamic_microgrids(out)
% Independent complete joint enumeration and fitted unilateral objectives.
root=fileparts(fileparts(mfilename('fullpath')));f=jsondecode(fileread(fullfile(root,'matlab/fixtures/microgrid_inputs.json')));a=rap.assets(f([f.n]==2));states=rap.product_states(3,4)+1;w=[.2 .6 .2];z=[-1 0 1];weights=prod(w(states),2);maxerr=zeros(1,5);count=0;
for fractions=[.1 .5 .9 .1;.1 .5 .9 .9]
 for step=[0 24 48 76 95]
  for scale=[0 .5 1 2]
   energy=a.capacity.*fractions';p=rap.physics(a,energy,step,scale);[mu,u,C,load]=rap.inputs(p,2.5);x=zeros(81,2);im=x;cost=x;ext=zeros(81,1);
   for j=1:81
    ix=states(j,:);x(j,:)=rap.draw(p,'export',ix);im(j,:)=rap.draw(p,'import',ix);cost(j,:)=5*(rap.draw(p,'charge',ix)+rap.draw(p,'discharge',ix));ext(j)=p.external_cond_mean(ix(1))+.02*scale*p.external_base*z(ix(4));
   end
   independent_load=weights'*(ext+sum(im,2));independent_mu=weights'*x;uc=zeros(1,2);
   for i=1:2,[v,ix]=sort(x(:,i));uc(i)=min(independent_mu(i),v(find(cumsum(weights(ix))>=.1,1)));end
   err=[max(abs(mu-independent_mu)),abs(load-independent_load),max(abs(u-uc)),0,0];[q,k]=rap.strategies(mu,u);
   for profile=[zeros(2,1) u' q']
    profile=profile';cash=settle(profile,0);pred=.25*(C+5*profile*(sum(mu)-sum(profile)));err(4)=max(err(4),max(abs(cash-pred)));
   end
   for j=1:3
    assert(max(abs(q(j,:)-rap.equilibrium(sum(mu),u,k(j))))<1e-9);
    for i=1:2
     if u(i)<1e-10,br=0;else
      vals=zeros(1,3);
      for h=1:3,trial=q(j,:);trial(i)=(h-1)*u(i)/2;[~,objective]=settle(trial,(k(j)-1)/10);vals(h)=objective(i);end
      A=2*(vals(3)-2*vals(2)+vals(1))/u(i)^2;B=(vals(3)-vals(1)-A*u(i)^2)/u(i);assert(A<0);br=min(u(i),max(0,-B/(2*A)));
     end
     err(5)=max(err(5),abs(br-q(j,i)));
    end
   end
   assert(max(weights'*(x<u-1e-12))<=.1+1e-12);assert(all(err<1e-9));maxerr=max(maxerr,err);count=count+1;
  end
 end
end
report=struct('status','passed','conditional_states',count,'joint_realizations',81*count,'max_errors',maxerr,'error_order',{{'export means','load mean','quantile caps','direct settlement','fitted best response'}});rap.save_report(out,'independent_microgrid_audit',report);
 function [cash,objective]=settle(q,beta)
  rt=20+5*(ext+sum(im,2)-sum(x,2));da=20+5*(independent_load-sum(q));cash=weights'*(.25*(da*q+rt.*(x-im-q)-cost));objective=cash-.25*beta*(weights'*(rt-da).^2);
 end
end
