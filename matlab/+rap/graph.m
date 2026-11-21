function g=graph(n,kind)
m=n+1;A=false(m);sh=[];
switch kind
 case 'star',A(1:n,m)=true;A(m,1:n)=true;
 case 'complete',A=~eye(m);
 otherwise
  if strcmp(kind,'ring'),sh=[1 m-1];else,s=max(2,floor(sqrt(m)));sh=unique(mod([1 -1 s -s],m));sh(sh==0)=[];end
  for i=1:m,A(i,mod(i-1+sh,m)+1)=true;end
end
[i,j]=find(triu(A,1));e=sortrows([i j]);d=sum(A,2);den=max(d)+1;
if strcmp(kind,'star'),num=n;rd=m;elseif strcmp(kind,'complete'),num=0;rd=1;else
 vals=zeros(1,m-1);for j=1:m-1,t=mod(j*sh,m);vals(j)=sum(min(t,m-t).^2);end
 rd=m*m*den;num=max(rd-8*min(vals),(2*max(d)-den)*m*m);
end
used=false(m,2*m);counts=[];
for j=1:size(e,1),c=find(~(used(e(j,1),:)|used(e(j,2),:)),1);used(e(j,:),c)=true;if numel(counts)<c,counts(c)=0;end;counts(c)=counts(c)+2;end
W=eye(m)-(diag(d)-double(A))/den;rho=max(abs(eig(W)-0)); %#ok<NASGU>
g=struct('kind',kind,'m',m,'edges',e,'shifts',sh,'denominator',den,'rho_num',num,'rho_den',rd,'rho_upper',num/rd,'colors',counts,'W',W);
end
