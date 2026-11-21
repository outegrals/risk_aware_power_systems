function k=binomial(N,p)
% Inverse CDF using base MATLAB betainc; vector remaining counts, no toolbox.
if p<=0,k=zeros(size(N));return;end
if p>=1,k=N;return;end
u=rand(size(N));lo=zeros(size(N));hi=N;
while any(lo<hi)
 mid=floor((lo+hi)/2);cdf=ones(size(N));ix=mid<N;cdf(ix)=betainc(1-p,N(ix)-mid(ix),mid(ix)+1);
 left=u<=cdf;hi(left)=mid(left);lo(~left)=mid(~left)+1;
end
k=lo;
end
