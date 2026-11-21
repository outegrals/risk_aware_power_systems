function y=mask(x,edges,masks,p)
y=x(:);
for e=1:size(edges,1),i=edges(e,1);j=edges(e,2);y(i)=y(i)+masks(e);y(j)=y(j)-masks(e);end
y=mod(y,p);
end
