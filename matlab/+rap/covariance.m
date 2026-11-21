function C=covariance()
s=[12 10 8 20];rho=.3*ones(4);rho(4,:)=.2;rho(:,4)=.2;rho(1:5:end)=1;C=(s'*s).*rho;
end
