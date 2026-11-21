function z=product_states(base,n)
c=cell(1,n);[c{:}]=ndgrid(0:base-1);z=zeros(base^n,n);
for i=1:n,z(:,i)=c{i}(:);end
end
