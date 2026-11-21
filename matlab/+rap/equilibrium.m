function q=equilibrium(X,u,k)
lo=0;hi=X;
for j=1:100
 r=(lo+hi)/2;
 if r+sum(min(u,k.*r))>X,hi=r;else,lo=r;end
end
q=min(u,k.*((lo+hi)/2));
end
